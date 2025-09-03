#ifndef CSVWRITER_HPP
#define CSVWRITER_HPP

/**
 * Version: August 2025
 * Author: Nicolás Loubet
 */

#include "ToolKit.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <iostream>
#include <filesystem>

class CSVWriter {
    std::ofstream file;
	std::string filename;
    char sep;

	private:
		// Generic conversion to string
		template<typename T>
		static std::string toString(const T& value) {
			std::ostringstream oss;
			oss << value;
			return oss.str();
		}

		// Backup existing file
		void backupExistingFile(const std::string& filename) {
			if(std::filesystem::exists(filename)) {
				std::cout << "Warning: File " << filename << " already exists. Backing up existing content." << std::endl;
				
				int backupNum= 1;
				std::string backupFilename;
				do {
					backupFilename= filename + "." + std::to_string(backupNum);
					backupNum++;
				} while(std::filesystem::exists(backupFilename));
				
				try {
					std::filesystem::rename(filename, backupFilename);
					std::cout << "Previous content backed up to " << backupFilename << std::endl;
				} catch(const std::filesystem::filesystem_error& e) {
					throw std::runtime_error("Error creating backup: " + std::string(e.what()));
				}
			}
		}

	public:
		CSVWriter(const std::string& filename, bool backup= true, char separator= ','): filename(filename), sep(separator) {
			if(!backup) return;
			backupExistingFile(filename);
			file.open(filename, std::ios::out | std::ios::trunc);
			if(!file.is_open()) throw std::runtime_error("Error: cannot open file " + filename);
		}

		~CSVWriter() {
			if(file.is_open()) file.close();
		}

		/**
		 * Writes the header of the CSV file
		 * @param headers vector of column names
		 */
		void writeHeader(const std::vector<std::string>& headers) {
			for(size_t i= 0; i < headers.size(); i++) {
				file << headers[i];
				if(i+1 < headers.size()) file << sep;
			}
			file << "\n";
		}

		/**
		 * Writes a row of the CSV file
		 * @param args list of column values
		 * @tparam Args types of the column values
		 */
		template<typename... Args>
		void writeRow(Args... args) {
			std::vector<std::string> values= {toString(args)...};
			for(size_t i= 0; i < values.size(); i++) {
				file << values[i];
				if(i+1 < values.size()) file << sep;
			}
			file << "\n";
		}

		/**
		 * Writes a row from a vector of generic type T
		 * @param row vector containing the column values
		 * @tparam T type of the elements in the vector
		 */
		template<typename T>
		void writeRow(const std::vector<T>& row) {
			for(size_t i= 0; i < row.size(); i++) {
				file << toString(row[i]);
				if(i+1 < row.size()) file << sep;
			}
			file << "\n";
		}

		/**
		 * Writes a normalized distribution in several columns
		 * @param bins vector of pointers to columns (each column has N_BINS values)
		 * @param N_BINS number of bins
		 * @param MIN_BINS minimum value
		 * @param MAX_BINS maximum value
		 * @param titles column titles (the first one is for the x axis)
		 * @param totals (optional) normalizers by column. If empty, they are calculated by summing bins.
		 */
		void writeDistribution(const std::vector<Real*>& bins, int N_BINS, Real MIN_BINS, Real MAX_BINS,
							   const std::vector<std::string>& titles, const std::vector<Real>& totals= {}) {
			if(titles.size() != bins.size()+1) throw std::runtime_error("writeDistribution: titles.size() must be bins.size()+1");

			std::vector<Real> norm(bins.size());
			if(totals.empty()) {
				for(size_t c= 0; c < bins.size(); c++)
					norm[c]= std::accumulate(bins[c], bins[c]+N_BINS, 0.0);
			} else {
				norm = totals;
			}

			writeHeader(titles);

			const float DX= (MAX_BINS - MIN_BINS) / N_BINS;

			for(int i= 0; i < N_BINS; i++) {
				file << (i*DX + MIN_BINS + DX/2.0);
				for(size_t c= 0; c < bins.size(); c++) {
					Real val= (norm[c] > 0 ? bins[c][i] / norm[c] : 0.0);
					file << sep << val;
				}
				file << "\n";
			}
		}

		/**
		 * Writes a generic table
		 * @param rows row names (first column)
		 * @param table pointer double with the data [i_row][i_col]
		 * @param nRows number of rows
		 * @param nCols number of columns
		 * @param colTitles column titles (must be nCols+1, including first)
		 */
		template<typename T>
		void writeTable(const std::vector<std::string>& rows, T** table, const int nRows, const int nCols, const std::vector<std::string>& colTitles) {
			if(rows.size() != nRows) throw std::runtime_error("writeTable: rows.size() != nRows");
			if(colTitles.size() != nCols+1) throw std::runtime_error("writeTable: colTitles.size() != nCols+1");

			writeHeader(colTitles);

			for(int i= 0; i < nRows; i++) {
				file << rows[i];
				for(int j= 0; j < nCols; j++)
					file << sep << table[i][j];
				file << "\n";
			}
		}

		/**
		 * Appends a single column to the right of an existing CSV file
		 * @param column vector containing the column values
		 * @param columnTitle title of the new column
		 * @param nRows expected number of rows
		 * @tparam T type of the elements in the column
		 */
		template<typename T>
		void appendColumn(const std::vector<T>& column, const std::string& columnTitle, size_t nRows) {
			if(column.size() != nRows) throw std::runtime_error("appendColumn: column.size() must equal nRows");

			std::ifstream inFile(filename, std::ios::in);
			if(!inFile.is_open()) throw std::runtime_error("Error: cannot read file " + filename);

			std::vector<std::string> lines;
			std::string line;
			while(std::getline(inFile, line))
				lines.push_back(line);
			inFile.close();

			if(lines.empty()) throw std::runtime_error("appendColumn: file is empty (no header to extend)");
			if(lines.size()-1 != nRows) throw std::runtime_error("appendColumn: nRows does not match number of data rows in file");
			if(file.is_open()) file.close();

			std::ofstream out(filename, std::ios::out | std::ios::trunc);
			if(!out.is_open()) throw std::runtime_error("Error: cannot open file " + filename);

			out << lines[0] << sep << columnTitle << "\n";
			for(size_t i= 1; i < lines.size(); i++)
				out << lines[i] << sep << toString(column[i-1]) << "\n";
			out.close();
		}

		/**
		 * Appends multiple columns to the right of an existing CSV file
		 * @param columns vector of vectors containing the column values
		 * @param columnTitles titles of the new columns
		 * @param nRows expected number of rows
		 * @tparam T type of the elements in the columns
		 */
		template<typename T>
		void appendColumns(std::vector<std::vector<T>> columns, const std::vector<std::string>& columnTitles, size_t nRows) {
			bool isColumns= (columns.size() == columnTitles.size());
			bool isRows= (columns.size() == nRows && !columns.empty() && columns[0].size() == columnTitles.size());

			if(isColumns) {
				for(const auto& column: columns) {
					if(column.size() != nRows) throw std::runtime_error("appendColumns: each column.size() must equal nRows");
				}
			} else if(isRows) {
				std::vector<std::vector<T>> transposed(columnTitles.size(), std::vector<T>(nRows));
				for(size_t i= 0; i < nRows; i++) {
					if(columns[i].size() != columnTitles.size()) throw std::runtime_error("appendColumns: row size mismatch with columnTitles");
					for(size_t j= 0; j < columnTitles.size(); j++) {
						transposed[j][i]= columns[i][j];
					}
				}
				columns= std::move(transposed);
			} else {
				throw std::runtime_error("appendColumns: dimensions do not match either (nCols x nRows) or (nRows x nCols)");
			}

			std::ifstream inFile(filename, std::ios::in);
			if(!inFile.is_open()) throw std::runtime_error("Error: cannot read file " + filename);

			std::vector<std::string> lines;
			std::string line;
			while(std::getline(inFile, line))
				lines.push_back(line);
			inFile.close();

			if(lines.empty()) throw std::runtime_error("appendColumns: file is empty (no header to extend)");
			if(lines.size()-1 != nRows) throw std::runtime_error("appendColumns: nRows does not match number of data rows in file");
			if(file.is_open()) file.close();

			std::ofstream out(filename, std::ios::out | std::ios::trunc);
			if(!out.is_open()) throw std::runtime_error("Error: cannot open file " + filename);

			out << lines[0];
			for(const auto& title: columnTitles)
				out << sep << title;
			out << "\n";

			for(size_t i= 1; i < lines.size(); i++) {
				out << lines[i];
				for(const auto& column: columns)
					out << sep << toString(column[i-1]);
				out << "\n";
			}

			out.close();
		}

};

#endif // CSVWRITER_HPP
