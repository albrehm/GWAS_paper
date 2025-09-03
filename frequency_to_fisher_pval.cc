#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <cmath>  

// Fisher's exact test from htslib
// https://github.com/samtools/htslib
// kfunc.h file from htslib
#include "kfunc.h"

using namespace std;

// Function to print detailed usage information
void printUsage(const char* progName) {
    cout << "Usage: " << progName << " [options] <csv_filename>\n"
         << "Options:\n"
         << "  -h, --help             Show this help message and exit\n"
         << "  -o, --output <file>    Specify an output file to save the processed data\n"
         << "\n"
         << "The CSV file should have exactly 8 columns per row, with the first line as a header.\n"
         << "The CSV file should be the output of grenedalf frequency, run with 2 input SAM/BAM.\n";
}

// Structure to represent a row of data (with 8 columns)
struct DataRow {
    string col1, col2, col3, col4, col5, col6, col7, col8;
};

bool importCSV(const string& filename, vector<DataRow>& data) {
    ifstream file(filename);  // Open the file

    // Check if the file was successfully opened
    if (!file.is_open()) {
        cerr << "Could not open the file: " << filename << endl;
        return false;
    }

    string line;
    bool isFirstLine = true;  // Flag to skip the header line

    while (getline(file, line)) {
        if (isFirstLine) {
            isFirstLine = false;  // Skip header
            continue;
        }

        stringstream ss(line);  // String stream to parse the line
        DataRow row;
        string temp;

        // Read each column into the DataRow struct
        getline(ss, row.col1, ',');
        getline(ss, row.col2, ',');
        getline(ss, row.col3, ',');
        getline(ss, row.col4, ',');
        getline(ss, row.col5, ',');
        getline(ss, row.col6, ',');
        getline(ss, row.col7, ',');
        getline(ss, row.col8, ',');

        data.push_back(row);  // Add the row to the data vector
    }

    file.close();  // Close the file
    return true;
}

// Function to print the data
void printData(const vector<DataRow>& data) {
    for (const auto& row : data) {
        cout << row.col1 << ", " << row.col2 << ", " << row.col3 << ", " 
             << row.col4 << ", " << row.col5 << ", " << row.col6 << ", " 
             << row.col7 << ", " << row.col8 << endl;
    }
}

double fisher_twosided_pval (int ref_1, int alt_1, int ref_2, int alt_2)  {
  double fisher_left_p, fisher_right_p, fisher_twosided_p;

  // from htslib 
  kt_fisher_exact(ref_1, alt_1, ref_2, alt_2,
    &fisher_left_p, &fisher_right_p, &fisher_twosided_p);

  return fisher_twosided_p;
}

// calculate a p-value for each row
void calculate_pvalues(const vector<DataRow>& data) {
    for (const auto& row : data) {

		  int ref_1 = std::stoi(row.col5);
		  int alt_1 = std::stoi(row.col6);
		  int ref_2 = std::stoi(row.col7);
		  int alt_2 = std::stoi(row.col8);

	     double pval = fisher_twosided_pval(ref_1, alt_1, ref_2, alt_2);
		  double log_pval = std::log10(pval);
		  double neg_log_pval = (log_pval == 0) ? 0 : -log_pval;
	     
		  // write results to stdout
		  // original input + -log10 of fisher p-val
        cout << row.col1 << ", " << row.col2 << ", " << row.col3 << ", " 
             << row.col4 << ", " << row.col5 << ", " << row.col6 << ", " 
             << row.col7 << ", " << row.col8 << ", " << neg_log_pval
				 << endl;
    }
}

int main(int argc, char* argv[]) {

    bool verbose = false;
    string outputFilename;
    bool listOnly = false;

    // read in command-line options
    int opt;
    while ((opt = getopt(argc, argv, "ho:")) != -1) {
        switch (opt) {
            case 'h':
                printUsage(argv[0]);
                return EXIT_SUCCESS;
				case 'o':
                outputFilename = optarg;  // Get the argument for -o option
                break;
            default:
                cerr << "Invalid option. Use -h for help.\n";
                printUsage(argv[0]);
                return EXIT_FAILURE;
        }
    }

	 // After processing options, we expect the CSV filename to be the last argument
    if (optind >= argc) {
        cerr << "Error: Missing required CSV filename argument.\n";
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

	 string filename = argv[optind];  // Get the CSV filename

    vector<DataRow> data;

    // Import the CSV data
    if (importCSV(filename, data)) {
        if (verbose) {
            cout << "Successfully imported CSV data from: " << filename << endl;
        }

	     } else {
        return EXIT_FAILURE;  // If CSV import failed
    }

	 // calculate pvalues and output to stdout
	 calculate_pvalues(data);  

    return EXIT_SUCCESS;
}


