// TOML reading - includes uninspected AI generated code
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "toml.hpp"
#include "../io.cpp"
// Include the io.cpp file to access the BinningConfig, RuntimeConfig, and AnalysisConfig structs




std::vector<double> load_binary_doubles(const std::string& filename, std::size_t expected_values) {
    /**
     * Load double-precision values from a binary file.
     *
     * @param filename:         Path to the binary file.
     * @param expected_values:  Expected number of values, or zero to skip validation.
     * @return:                 The values read from the file.
     * @throws:                 std::runtime_error If the file is invalid, cannot be read, or has an
     *                          unexpected number of values.
     */
    
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Could not open binary file: " + filename);
    }

    const std::streamsize file_size = file.tellg();
    if (file_size < 0 || file_size % static_cast<std::streamsize>(sizeof(double)) != 0) {
        throw std::runtime_error("Binary file does not contain a whole number of doubles: " + filename);
    }

    const std::size_t value_count = static_cast<std::size_t>(file_size) / sizeof(double);
    if (expected_values != 0 && value_count != expected_values) {
        throw std::runtime_error("Unexpected number of values in " + filename +
                                 ": expected " + std::to_string(expected_values) +
                                 ", found " + std::to_string(value_count));
    }

    std::vector<double> values(value_count);
    file.seekg(0, std::ios::beg);
    file.read(reinterpret_cast<char*>(values.data()), file_size);
    if (!file) {
        throw std::runtime_error("Failed while reading binary file: " + filename);
    }
    return values;
}

toml::table read_toml_file(const std::string& filename) {
    /**
    * Read a TOML file and return its contents as a toml::table.
    *
    * @param filename:         Path to the TOML file.
    * @return:                 toml::table containing the parsed TOML data.
    * @throws:                 std::runtime_error If the file cannot be opened or parsed.
    */
    try {
        toml::table config = toml::parse_file(filename);
        return config;
    } catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        throw; // Rethrow the exception to be handled by the caller
    }
}

int calc_N_components(const toml::array& fields){
    /**
    * Calculate the total number of components based on the provided fields.
    *
    * @param fields:           toml::array containing the field definitions.
    * @return:                 Total number of components.
    */
    int N_components = 3; // start with 3 total components for the position field (x, y, z)

    int i = 0;

    // loop through each field in the fields array
    for (const toml::node& field_node : fields) {
        // Access the field as a table
        const toml::table* field = field_node.as_table();
        // Check if the field is valid and contains the "components" key
        if (!field) {
            std::cerr << "Warning: Field " << i << " is not a table. Assuming 1 component.\n";
            N_components += 1;
        } else if (field->contains("components")) {
            // Get the "components" value, defaulting to 1 if not specified
            const auto* components_node = field->get_as<int64_t>("components");
            const int64_t components = components_node ? components_node->get() : 1;
            
            // Add the number of components for this field to the total
            N_components += static_cast<int>(components);
        } else {
            std::cerr << "Warning: Field " << i << " does not specify 'components'. Assuming 1 component.\n";
            N_components += 1; // Default to 1 if not specified
        }
        ++i;
    }

    return N_components;
}

int calc_N_histograms(const toml::array& fields){
    /**
    * Calculate the total number of histograms based on the provided fields.
    *
    * @param fields:           toml::array containing the field definitions.
    * @return:                 Total number of histograms.
    */
    int N_histograms = 0; // start with 0 histograms

    int i = 0;

    // loop through each field in the fields array
    for (const toml::node& field_node : fields) {
        // Access the field as a table
        const toml::table* field = field_node.as_table();

        // Get the "statistics" and "orders" arrays, defaulting to nullptr if not specified
        const toml::array* statistics = field ? field->get_as<toml::array>("statistics") : nullptr;
        const toml::array* orders = field ? field->get_as<toml::array>("orders") : nullptr;

        // Calculate the number of histograms for this field based on the sizes of the "statistics" and "orders" arrays
        if (statistics && orders) {
            N_histograms += static_cast<int>(statistics->size() * orders->size());
        } 
        
        // handle lack of specification for either "statistics" or "orders"
        else {
            int stats_count;
            int orders_count;

            // If either "statistics" or "orders" is not specified, assume 1 for the missing one
            if (!statistics){
                std::cerr << "Warning: Field " << i << " does not specify 'statistics'. Assuming scalar and thus 1 statistic.\n";
                stats_count = 1;
            }
            else{
                stats_count = static_cast<int>(statistics->size());
            }

            if (!orders){
                std::cerr << "Warning: Field " << i << " does not specify 'orders'. Assuming 1 order and thus 1 histogram per statistic.\n";
                orders_count = 1;
            }
            else{
                orders_count = static_cast<int>(orders->size());
            }

            std::cerr << "Warning: Field " << i << " (Property = " << field->get_as<std::string>("name")->get() << ") does not specify 'statistics' and/or 'orders'.\n";
            std::cerr << "--> Assuming "<< stats_count <<" statistic and "<< orders_count <<" orders,";
            std::cerr << "thus adding "<< stats_count * orders_count <<" histograms for this field.\n";
            N_histograms += stats_count * orders_count;
        }
        ++i;
    }
    return N_histograms;
}



int main(int argc, char* argv[]) {

    // Start by loading the configuration from the TOML file specified in the command line arguments
    toml::table config;
    try {
        config = read_toml_file(argv[1]);
        std::cout << config << "\n";
    }
    catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        return 1;
    }

    // Print the "input" section of the configuration to the console for verification
    std::cout << config["input"] << "\n-------\n";

    // From the config file, treat provided "fields" as an array of tables
    const toml::array* fields = config["field"].as_array();

    // Check it is a table
    if (!fields) {
        std::cerr << "Expected 'field' to be an array of tables.\n";
        return 1;
    }

    // Load the positions from the binary file specified in the configuration
    const auto positions_filename = config["input"]["positions"].value<std::string>();
    if (!positions_filename) {
        std::cerr << "Expected input.positions to be a string.\n";
        return 1;
    }

    const std::vector<double> position_values = load_binary_doubles(*positions_filename, 0);
    
    // Ensure 3D position has been provided for each particle
    // TODO: Add array handling to support 2D or 1D positions in the future - this will add zeros in place of missing dimensions as long as 1 is provided
    if (position_values.size() % 3 != 0) {
        std::cerr << "The positions file must contain 3 doubles per particle.\n";
        return 1;
    }

    // Calculate the number of particles, components, and histograms based on the loaded data and configuration
    /// -- Components counts e.g. x,y,z,vx,vy,vz, etc. separately
    /// -- Histograms track counts and statistic sums for each output required
    const std::size_t N_part = position_values.size() / 3;
    const int N_components = calc_N_components(*fields);
    const int N_histograms = calc_N_histograms(*fields);

    
    
    // create a vector of vectors to hold the input data
    /// Row-major input data: columns 0-2 are positions, followed by configured fields.
    std::vector<std::vector<double>> input_data(N_part,  // N_part rows
                                                std::vector<double>(static_cast<std::size_t>(N_components), // N_components columns
                                                0.0)); // initialize all values to 0.0
    /// Fill the first three columns with the position data
    for (std::size_t particle = 0; particle < N_part; ++particle) {
        for (std::size_t component = 0; component < 3; ++component) {
            input_data[particle][component] = position_values[particle * 3 + component];
        }
    }

    // create a vector of vectors to hold the output histograms
    //TODO: Add ability to set N_bins from the config file, currently hardcoded to 200
    std::vector<std::vector<double>> output_histograms(200, // 200 histogram bins
                                                       std::vector<double>(static_cast<std::size_t>(N_histograms), // N_histograms columns
                                                       0.0)); // initialize all values to 0.0

    // Prepare to populate arrays

    /// To ease looping over properties, initialize the component offset to 3,
    /// since the first three components are for position (x, y, z)
    std::size_t component_offset = 3;

    /// Get the number of specified fields using fields->size() and print it to the console
    std::cout << "Number of fields: " << fields->size() << "\n";
    
    // Populate the "input_data" array with the values from each field specified in the configuration

    /// For each entry in fields
    for (const toml::node& field_node : *fields) {
        
        // access it as a table
        const toml::table* field = field_node.as_table();

        // error if not a table
        if (!field) {
            std::cerr << "Expected every field entry to be a table.\n";
            return 1;
        }

        /// pull out each defined property for the property
        // -- property name
        const auto name = field->get_as<std::string>("name");  
        // -- binary file containing the property values
        const auto filename = field->get_as<std::string>("file");  
        // -- number of components for the property (e.g., 1 for scalar, 3 for vector)
        const auto components = field->get_as<int64_t>("components");  
        // -- array of statistics to compute for the property (e.g., VSF, VSF perp, VSF parallel)
        const toml::array* statistics = field->get_as<toml::array>("statistics"); 
        // -- array of orders to compute for the property (e.g., 2nd order) 
        const toml::array* orders = field->get_as<toml::array>("orders");

        // error if any of the required properties are missing
        if (!name || !filename || !components) {
            std::cerr << "Each field must contain name, file, and components.\n";
            return 1;
        }

        // error if the number of components is less than 1
        const int64_t field_components = components->get();
        if (field_components < 1) {
            std::cerr << "Field " << name->get() << " must have at least one component.\n";
            return 1;
        }

        // Load the binary data for the field and populate the corresponding columns in the input_data array
        const std::vector<double> field_values = load_binary_doubles(filename->get(), // filename pulled from toml config
                                                                     N_part * static_cast<std::size_t>(field_components)); // expected number of values based on the number of particles and components

        // Populate the input_data array with the loaded field values
        for (std::size_t particle = 0; particle < N_part; ++particle) {
            for (std::size_t component = 0;  // For component
                 component < static_cast<std::size_t>(field_components); // in range N_components for this field
                 ++component) { // incrementing by 1
                // Set value for particle j and component k in the input_data array
                input_data[particle][component_offset + component] =
                    field_values[particle * static_cast<std::size_t>(field_components) + component];
            }
        }
        // increase the component offset by the number of components already allocated to input array
        component_offset += static_cast<std::size_t>(field_components);

        // Print details about the field to the console for verification
        std::cout << "Field details:\n";
        std::cout << "Name: " << *name << "\n";
        std::cout << "File: " << *filename << "\n";
        std::cout << "Components: " << field_components << "\n";
        std::cout << "Statistics:";
        if (statistics) {
            for (const toml::node& statistic : *statistics) {
                if (const auto value = statistic.value<std::string>()) {
                    std::cout << " " << *value;
                }
            }
        }
        std::cout << "\nOrders:";
        if (orders) {
            for (const toml::node& order : *orders) {
                if (const auto value = order.value<int64_t>()) {
                    std::cout << " " << *value;
                }
            }
        }
        std::cout << "\n";


    }

    // Verify that the total number of components matches the expected number based on the configuration
    if (component_offset != static_cast<std::size_t>(N_components)) {
        std::cerr << "Loaded field components do not match N_components.\n";
        return 1;
    }
    std::cout << "Total number of components: " << N_components << "\n";

    std::cout << "Total number of histograms: " << N_histograms << "\n";
    std::cout << "Loaded input array: " << N_part << " x " << N_components << "\n";
    std::cout << "Allocated output array: 200 x " << N_histograms << "\n";

    

    return 0;
}