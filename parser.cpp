#include "parser.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

std::vector<double> readCSV(std::string filename, int index,
                            bool printing) {

    // Setup file
    std::ifstream myfile;
    myfile.open(filename);
    std::string line;
    std::getline(myfile, line);

    // Setup output and temps
    std::vector<double> output;
    std::vector<std::string> temps;
    std::string temp;

    // Parsing
    if (myfile.is_open()) {
        while (myfile) {
            std::getline(myfile, line);
            if (line.length() == 0) {
                continue;
            }
            std::stringstream ss(line);
            while (std::getline(ss, temp, ',')) {
                temps.push_back(temp);
            }
            output.push_back(std::stod(temps[index]));
            temps = std::vector<std::string>();
        }
    }
    if (printing) {
        std::cout << "Finish parsing!" << std::endl;
        std::cout << "Parsed vector of length " << output.size() << std::endl;
    }
    return output;
}

void saveJSON(std::vector<double> input, std::string filename,
              std::string plot_name, std::string x_label ,
              std::string y_label, std::string x_type,
              std::string x_min,
              std::string x_max) {
    //         f = fopen(filename.c_str(), "w+");
    // fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\"
    // height = \"1000\">\n"); fprintf(f, "<g>\n");
    FILE *f = fopen(filename.c_str(), "w");
    fprintf(f, "{\n\t");
    fprintf(f, "\"Data\": [");
    for (int i = 0; i < input.size(); i++) {
        fprintf(f, "%3.5f", input[i]);
        if (i < input.size() - 1) {
            fprintf(f, ", ");
        }
    }
    fprintf(f, "],\n\"Size\": %zu,\n", input.size());
    const std::string name = "\"Plot name\": \"" + plot_name + "\",\n";
    fprintf(f, "%s", name.c_str());
    fprintf(f, "%s", ("\"X type\": \"" + x_type + "\",\n").c_str());
    fprintf(f, "%s", ("\"X min\": \"" + x_min + "\",\n").c_str());
    fprintf(f, "%s", ("\"X max\": \"" + x_max + "\",\n").c_str());
    fprintf(f, "%s", ("\"X label\":\"" + x_label + "\",\n").c_str());
    fprintf(f, "%s", ("\"Y label\": \"" + y_label + "\"\n").c_str());
    fprintf(f, "}\n\t");
}

void plotData(std::vector<double> input, std::string plotname,
              std::string filename, std::string plot_name,
              std::string x_label, std::string y_label,
              std::string x_type, std::string x_min,
              std::string x_max) {
    saveJSON(input, filename, plot_name, x_label, y_label, x_type, x_min,
             x_max);
    system(("python plotting.py " + filename + " " + plotname).c_str());
}

int test_parser() {
    std::string filename = "DailyDelhiClimateTrain.csv";
    std::vector<double> data = readCSV(filename, 1);
    for (int i = 0; i < data.size(); i++) {
        std::cout << data[i] << std::endl;
    }
    plotData(data);
    return 0;
}
