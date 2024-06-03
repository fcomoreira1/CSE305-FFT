#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

std::vector<double> readCSV(std::string filename, int index, bool printing = false) {
    
    // Setup file
    std::ifstream myfile;
    myfile.open(filename);
    std::string line;
    std::getline (myfile, line);

    // Setup output and temps
    std::vector<double> output;
    std::vector<std::string> temps;
    std::string temp;

    // Parsing
    if (myfile.is_open()) {
        while (myfile) {
            std::getline (myfile, line);
            if (line.length() == 0) {continue;}
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

void saveJSON(std::vector<double> input, std::string filename = "data.json", std::string plot_name = "", std::string x_label = "Date", std::string y_label = "Temperature", 
              std::string x_type = "date", std::string x_min = "2013-01-01", std::string x_max = "2017-01-01") {
        //         f = fopen(filename.c_str(), "w+");
        // fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        // fprintf(f, "<g>\n");
    FILE* f = fopen(filename.c_str(), "w+");
    fprintf(f, "{\n\t");
    fprintf(f, "Data: [");
    for (int i = 0; i < input.size(); i ++) {
        fprintf(f, "%3.5f, ", input[i]);
    }
    fprintf(f, "],\nSize: %d,\n", input.size());
    const std::string name = "Plot name: " + plot_name + ",\n";
    fprintf(f, name.c_str());
    fprintf(f, ("X type: "+x_type+",\n").c_str());
    fprintf(f, ("X min: "+x_min+",\n").c_str());
    fprintf(f, ("X max: "+x_max+",\n").c_str());
    fprintf(f, ("X label: "+x_label+",\n").c_str());
    fprintf(f, ("Y label: "+x_label+",\n").c_str());
}

void plotData(std::vector<double> input, std::string plotname = "plot.png", std::string filename = "data.json", std::string plot_name = "", std::string x_label = "Date", std::string y_label = "Temperature", 
              std::string x_type = "date", std::string x_min = "2013-01-01", std::string x_max = "2017-01-01") {
    saveJSON(input, filename, plot_name, x_label, y_label, x_type, x_min, x_max);
    system(("python plotting.py "+filename+" "+plotname).c_str());
}

int test() {
    std::string filename = "DailyDelhiClimateTrain.csv";
    std::vector<double> data = readCSV(filename, 1);
    for (int i = 0; i < data.size(); i ++) {
        std::cout << data[i] << std::endl;
    }
    plotData(data);
    return 0;
}