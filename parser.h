#pragma once
#include <string>
#include <vector>
std::vector<double> readCSV(std::string filename, int index,
                            bool printing = false);
void saveJSON(std::vector<double> input, std::string filename = "data.json",
              std::string plot_name = "plot", std::string x_label = "Date",
              std::string y_label = "Temperature", std::string x_type = "date",
              std::string x_min = "2013-01-01",
              std::string x_max = "2017-01-01");
void saveFile(std::vector<double> input, std::string filename = "data.txt",
              std::string plot_name = "plot", std::string x_label = "Date",
              std::string y_label = "Temperature", std::string x_type = "date",
              std::string x_min = "2013-01-01",
              std::string x_max = "2017-01-01");
void plotData(std::vector<double> input, std::string plotname = "plot.png",
              std::string filename = "data.json", std::string plot_name = "plot",
              std::string x_label = "Date", std::string y_label = "Temperature",
              std::string x_type = "date", std::string x_min = "2013-01-01",
              std::string x_max = "2017-01-01");
int test_parser();
