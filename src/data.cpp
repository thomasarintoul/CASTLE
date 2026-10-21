#include <algorithm>  // std::find_if
#include <iterator>   // std::distance
#include <cstddef>    // std::size_t


enum class Statistic {
    ScalarDifference,   // Used for 1D scalar properties (e.g., temperature, density, |B|)
    // Note: When the python frontend is implemented, we will need to enable it to create e.g. |B| if B is provided for ScalarDifference
    VectorDifference,   // Used for vector properties (e.g., velocity, B field)
    Parallel,           // Used for vector properties, computes the component of the difference along the separation vector
    Perpendicular       // Used for vector properties, computes the component of the difference perpendicular to the separation vector
};

enum class PropertyType {
    Scalar, // 1 component
    Vector  // 3 components
};

struct InputProperty {
    std::string name;               // Name of the property (e.g., "velocity", "temperature")
    PropertyType type;              // Type of the property (e.g., "scalar", "vector")
    std::vector<std::size_t> components;                 // Component "columns" in the input array (1 for scalar, 3 for vector)
};

struct OutputChannel{
    std::size_t input_property_index;  // Index of the input property in the input array
    Statistic statistic;                // Statistic to compute for this property (e.g., ScalarDifference, VectorDifference, Parallel, Perpendicular)
    int order;                          // Order of the statistic (e.g., 2 for second order, 3 for third order)
    std::size_t output_histogram_index; // Index of the output histogram "column" in the output array
};

void Properties(std::vector<InputProperty>& properties, const toml::array& fields) {
    // Input: vector of InputProperty structs already containing position;
    // and the toml::array of fields from the configuration

    std::size_t component_offset = 3; // Start after the first three components for position (x, y, z)

    // Loop through each field in the fields array and create an InputProperty struct for each
    for (const toml::node& field_node : fields) {
        // Each field must be represented by a TOML table
        const toml::table* field = field_node.as_table();
        if (!field) {
            throw std::runtime_error("Expected every field entry to be a table.");
        }

        // Populate the InputProperty struct with the relevant information from the TOML configuration
        const auto* name = field->get_as<std::string>("name");
        const auto* components = field->get_as<std::int64_t>("components");

        if (!name || !components || components->get() < 1) {
            throw std::runtime_error("Invalid or missing field name or component count.");
        }

        InputProperty property; // Create a new InputProperty struct for this field
        property.name = name->get();  // Set the property name
        property.type = (components->get() == 1) ? PropertyType::Scalar : PropertyType::Vector;  // Determine the property type based on the number of components

        // Based on the number of components specified, assign the appropriate component
        // indices for the input array
        for (std::size_t i = 0; i < static_cast<std::size_t>(components->get()); ++i) {
            property.components.push_back(component_offset + i);
        }

        // Then update the component offset for the next field
        component_offset += static_cast<std::size_t>(components->get());

        // Add the populated InputProperty struct to the properties vector
        properties.push_back(property);
    }
}

std::vector<OutputChannel> map_output_channels(const std::vector<InputProperty>& properties, const toml::array& fields) {
    // Create an array in which to store output channels
    // One output channel per histogram which corresponds to a unique property-statistic-order combination
    std::vector<OutputChannel> output_channels;

    // We assign histogram columns consecutively as channels are created
    std::size_t output_histogram_index = 0;

    // Loop over each field in the array
    for (const toml::node& field_node : fields) {

        // Each field must be represented by a TOML table
        const toml::table* field = field_node.as_table();
        
        // If the field is not a table, error and move on
        if (!field) {
            std::cerr << "Expected every field entry to be a table.\n";
            continue; // Skip this field and continue with the next
        }

        // Read the field name so it can be matched to an input property
        const auto name = field->get_as<std::string>("name");
        
        if (!name) {
            std::cerr << "Field does not specify a name.\n";
            continue; // Skip this field and continue with the next
        }

        // Find the corresponding input property by its configured name
        auto input_index = std::find_if(properties.begin(), properties.end(), [&](const InputProperty& prop) {
            return prop.name == *name;
        });
        
        // If we get to the end of the properties vector with no match, error and move on
        if (input_index == properties.end()) {
            std::cerr << "No matching InputProperty found for field: " << *name << "\n";
            continue; // Skip this field and continue with the next
        }

        // Get index in input vector from the iterator
        std::size_t input_property_index = std::distance(properties.begin(), input_index);

        // Retrieve statistics and orders from the TOML field.
        const toml::array* statistics = field->get_as<toml::array>("statistics");
        const toml::array* orders = field->get_as<toml::array>("orders");

        if (!orders) {
            throw std::runtime_error(
                "Field " + name->get() + " does not specify 'orders'."
            );
        }

        // Provide a local default for scalar properties.
        toml::array default_statistics;

        if (!statistics) {
            if (properties[input_property_index].type == PropertyType::Scalar) {
                std::cerr << "Field " << name->get()
                        << " does not specify 'statistics'. "
                        << "Assuming ScalarDifference.\n";
                        default_statistics.push_back("ScalarDifference");
                        statistics = &default_statistics;
            } else {
                std::cerr << "Field " << name->get()
                        << " does not specify 'statistics'." << "\n";
                continue; // Skip this field and continue with the next
            }
        }
        
    

        // We create individual channels for every statistic
        for (const toml::node& statistic_node : *statistics) {
            Statistic statistic; // variable to hold the statistic type

            // Convert the statistic name from TOML into the corresponding enum value.
            if (const auto stat_str = statistic_node.value<std::string>()) {
                if (*stat_str == "ScalarDifference" || *stat_str == "scalar") {
                    statistic = Statistic::ScalarDifference;
                } else if (*stat_str == "VectorDifference" || *stat_str == "vector") {
                    statistic = Statistic::VectorDifference;
                } else if (*stat_str == "Parallel" || *stat_str == "par" || *stat_str == "parallel" || *stat_str == "para") {
                    statistic = Statistic::Parallel;
                } else if (*stat_str == "Perpendicular" || *stat_str == "perp" || *stat_str == "perpendicular") {
                    statistic = Statistic::Perpendicular;
                } else {
                    std::cerr << "Unknown statistic: " << *stat_str << "\n";
                    continue; // Skip this statistic and continue with the next
                }
            } else {
                std::cerr << "Statistic is not a string for field: " << *name << "\n";
                continue; // Skip this statistic and continue with the next
            }

            // Within each statistic, we create individual channels for every requested order
            for (const toml::node& order_node : *orders) {
                if (const auto order_value = order_node.value<int64_t>()) {
                    // A valid statistic/order pair maps to one output histogram channel.
                    OutputChannel channel;  // Create new OutputChannel struct for this property-statistic-order combination

                    // populate the channel with the relevant information
                    channel.input_property_index = input_property_index;
                    channel.statistic = statistic;
                    channel.order = static_cast<int>(*order_value);
                    // Reserve the next available output column for this channel.
                    channel.output_histogram_index = output_histogram_index++;
                    
                    // Add the populated OutputChannel struct to the output_channels vector
                    output_channels.push_back(channel);
                } else {
                    std::cerr << "Order is not an integer for field: " << *name << "\n";
                }
            }
        }
    }
    return output_channels;
}