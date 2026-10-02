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
    try {
        toml::table config = toml::parse_file(filename);
        return config;
    } catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        throw; // Rethrow the exception to be handled by the caller
    }
}

int calc_N_components(const toml::array& fields){
    int N_components = 3; // start with 3 total components for the position field (x, y, z)

    int i = 0;
    for (const toml::node& field_node : fields) {
        const toml::table* field = field_node.as_table();
        if (!field) {
            std::cerr << "Warning: Field " << i << " is not a table. Assuming 1 component.\n";
            N_components += 1;
        } else if (field->contains("components")) {
            const auto* components_node = field->get_as<int64_t>("components");
            const int64_t components = components_node ? components_node->get() : 1;
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
    int N_histograms = 0; // start with 0 histograms

    int i = 0;
    for (const toml::node& field_node : fields) {
        const toml::table* field = field_node.as_table();

        const toml::array* statistics = field ? field->get_as<toml::array>("statistics") : nullptr;
        const toml::array* orders = field ? field->get_as<toml::array>("orders") : nullptr;
        if (statistics && orders) {
            N_histograms += static_cast<int>(statistics->size() * orders->size());
        } else {
            int stats_count;
            int orders_count;

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
    toml::table config;
    try {
        config = read_toml_file(argv[1]);
        std::cout << config << "\n";
    }
    catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        return 1;
    }

    std::cout << config["input"] << "\n-------\n";

    // treat fields as an array of tables
    const toml::array* fields = config["field"].as_array();
    if (!fields) {
        std::cerr << "Expected 'field' to be an array of tables.\n";
        return 1;
    }

    const auto positions_filename = config["input"]["positions"].value<std::string>();
    if (!positions_filename) {
        std::cerr << "Expected input.positions to be a string.\n";
        return 1;
    }

    const std::vector<double> position_values = load_binary_doubles(*positions_filename, 0);
    if (position_values.size() % 3 != 0) {
        std::cerr << "The positions file must contain 3 doubles per particle.\n";
        return 1;
    }

    const std::size_t N_part = position_values.size() / 3;
    const int N_components = calc_N_components(*fields);
    const int N_histograms = calc_N_histograms(*fields);

    // Row-major input data: columns 0-2 are positions, followed by configured fields.
    std::vector<std::vector<double>> input_data(
        N_part, std::vector<double>(static_cast<std::size_t>(N_components), 0.0));
    for (std::size_t particle = 0; particle < N_part; ++particle) {
        for (std::size_t component = 0; component < 3; ++component) {
            input_data[particle][component] = position_values[particle * 3 + component];
        }
    }

    // 200 bins/rows by one column per configured histogram.
    std::vector<std::vector<double>> output_histograms(
        200, std::vector<double>(static_cast<std::size_t>(N_histograms), 0.0));

    std::size_t component_offset = 3;

    // Get the number of specified fields using fields->size() and print it to the console
    std::cout << "Number of fields: " << fields->size() << "\n";
    
    // for each entry in fields
    for (const toml::node& field_node : *fields) {
        
        // access it as a table
        const toml::table* field = field_node.as_table();

        if (!field) {
            std::cerr << "Expected every field entry to be a table.\n";
            return 1;
        }

        // pull out each defined property for the property
        const auto name = field->get_as<std::string>("name");
        const auto filename = field->get_as<std::string>("file");
        const auto components = field->get_as<int64_t>("components");
        const toml::array* statistics = field->get_as<toml::array>("statistics");
        const toml::array* orders = field->get_as<toml::array>("orders");

        if (!name || !filename || !components) {
            std::cerr << "Each field must contain name, file, and components.\n";
            return 1;
        }

        const int64_t field_components = components->get();
        if (field_components < 1) {
            std::cerr << "Field " << name->get() << " must have at least one component.\n";
            return 1;
        }

        const std::vector<double> field_values = load_binary_doubles(
            filename->get(), N_part * static_cast<std::size_t>(field_components));
        for (std::size_t particle = 0; particle < N_part; ++particle) {
            for (std::size_t component = 0;
                 component < static_cast<std::size_t>(field_components); ++component) {
                input_data[particle][component_offset + component] =
                    field_values[particle * static_cast<std::size_t>(field_components) + component];
            }
        }
        component_offset += static_cast<std::size_t>(field_components);

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