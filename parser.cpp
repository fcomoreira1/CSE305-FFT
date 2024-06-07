#include "parser.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

std::vector<double> readCSV(std::string filename, int index, bool printing) {

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
              std::string plot_name, std::string x_label, std::string y_label,
              std::string x_type, std::string x_min, std::string x_max) {
    //         f = fopen(filename.c_str(), "w+");
    // fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\"
    // height = \"1000\">\n"); fprintf(f, "<g>\n");
    FILE *f = fopen(filename.c_str(), "w");
    fprintf(f, "{\n\t");
    fprintf(f, R"("Data": [)");
    for (int i = 0; i < input.size(); i++) {
        fprintf(f, "%3.5f", input[i]);
        if (i < input.size() - 1) {
            fprintf(f, ", ");
        }
    }
    fprintf(f, "],\n\t\"Size\": %zu,\n\t", input.size());
    const std::string name = "\"Plot name\": \"" + plot_name + "\",\n\t";
    fprintf(f, "%s", name.c_str());
    fprintf(f, "%s", ("\"X type\": \"" + x_type + "\",\n\t").c_str());
    fprintf(f, "%s", ("\"X min\": \"" + x_min + "\",\n\t").c_str());
    fprintf(f, "%s", ("\"X max\": \"" + x_max + "\",\n\t").c_str());
    fprintf(f, "%s", ("\"X label\": \"" + x_label + "\",\n\t").c_str());
    fprintf(f, "%s", ("\"Y label\": \"" + y_label + "\"\n").c_str());
    fprintf(f, "}");
}

void saveFile(std::vector<double> input, std::string filename,
              std::string plot_name, std::string x_label, std::string y_label,
              std::string x_type, std::string x_min, std::string x_max) {
    FILE *f = fopen(filename.c_str(), "w");
    for (int i = 0; i < input.size(); i++) {
        fprintf(f, "%3.5f\n", input[i]);
    }
    fprintf(f, "\n%zu,%s,%s,%s,%s,%s,%s", input.size(), plot_name.c_str(),
            x_type.c_str(), x_min.c_str(), x_max.c_str(), x_label.c_str(),
            y_label.c_str());
}

void plotData(std::vector<double> input, std::string plotname,
              std::string filename, std::string plot_name, std::string x_label,
              std::string y_label, std::string x_type, std::string x_min,
              std::string x_max) {
    std::string args = "";
    for (int i = 0; i < input.size(); i++) {
        args.append(std::to_string(input[i]));
        args.append(" ");
    }
    args.append(std::to_string(input.size()));
    args.append(" ");
    args.append(plot_name);
    args.append(" ");
    args.append(x_type);
    args.append(" ");
    args.append(x_min);
    args.append(" ");
    args.append(x_max);
    args.append(" ");
    args.append(x_label);
    args.append(" ");
    args.append(y_label);
    args.append(" ");
    args.append(plotname);
    args.append(" ");

    // FILE *f = fopen(filename.c_str(), "w");
    // fprintf(f, args.c_str());
    // std::cout << "Finished generating file!" << std::endl;

    std::ofstream fs;
    fs.open("data/_in_data.txt",
            std::ofstream::in | std::ofstream::out | std::ofstream::trunc);
    fs << args;
    fs.close();
    std::cout << "Finished generating file!" << std::endl;
}

int test_parser() {
    std::string filename = "data/DailyDelhiClimateTrain.csv";
    std::vector<double> data = readCSV(filename, 1);
    // saveFile(data, "data.txt");
    // system("python plotting.py data.txt");
    plotData(data, "data/plot.png", "data/data.txt");
    system("python plotting.py data/_in_data.txt");
    return 0;
}
