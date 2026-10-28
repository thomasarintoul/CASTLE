#include "data.hpp"
#include "io.hpp"
#include "config.hpp"

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
        property.filename = field->get_as<std::string>("file")->get(); // Set the filename for the property
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

void populate_input_data(std::vector<double>& input_data, const std::vector<InputProperty>& properties, std::size_t N_part){
    // set component offset to not overwrite position data
    int component_offset = 3; // Start after the first three components for position (x, y, z)

    // Loop through each prop in the properties vector (skipping position)
    for (int i = 0; i < properties.size(); ++i) {
        // Skip the first property (position) since it's already populated
        if (properties[i].name == "position") {
            continue;
        }

        // Print the property being processed
        std::cout << "Populating input_data for property: " << properties[i].name << "\n";

        // Read data using load_binary_doubles
        const std::vector<double> values = load_binary_doubles(properties[i].filename, 
                                                               N_part * properties[i].components.size());

        // if size / N_part != components.size(), error
        if (values.size() / N_part != properties[i].components.size()) {
            throw std::runtime_error("Mismatch between number of values in file " + properties[i].filename +
                                     " and expected number of components for property " + properties[i].name);
        }
        else{
            // if size / N_part == components.size(), populate input_data with the values in the correct columns
            // Do this by looping over components to populate the correct columns in input_data for each particle
            // (this should work for any number of components)
            for (std::size_t component : properties[i].components){ // for each component listed in properties[i]
                std::cout << "Populating component " << component << " with property " << properties[i].name << "\n";
                std::cout << "Component offset: " << component_offset << "\n";
                std::cout << "Using column " << component - component_offset << " from file " << properties[i].filename << "\n";
                std::cout << "-------" << "\n";

                for (std::size_t particle = 0; particle < N_part; ++particle) { // loop over particles
                    // for ith particle in jth component, values is at value[N_part * (j - N_allocated_components) + i]
                    input_data[component * N_part + particle] = values[N_part * (component - component_offset) + particle];
                }
            }
            component_offset += properties[i].components.size(); // update component offset for next property
            }

    }    


}
