//
// Created by Boyuan Shi on 18/09/2025.
//

#ifndef FILE_COMBINATION_H
#define FILE_COMBINATION_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <type_traits>
#include <cctype>
#include <future>
#include <iomanip>
#include <future>
#include <filesystem>
#include <regex>
#include <map>
#include <set>
#include <thread>
#include <array>
#include <type_traits>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <cctype>
#include <future>
#include <iomanip>
#include <filesystem>
#include <regex>
#include <map>
#include <set>
#include <thread>
#include <array>
#include "parameters.h"

inline void combine_files_RG_MCMC(int seconds) {

    // Combine output files by order number
    std::this_thread::sleep_for(std::chrono::seconds(seconds));
    std::cout << "Starting file combination..." << std::endl;

    // Map to store files by order number with their job indices
    std::map<int, std::vector<std::pair<std::string, int>>> regular_files_by_order;  // path, job_index
    std::map<int, std::vector<std::pair<std::string, int>>> block_std_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> seed_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> sigma_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> accept_ratio_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> direction_accept_ratio_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> tracked_arg_files_by_order;
    std::map<int, std::vector<std::tuple<std::string, int, int>>> histogram_files_by_order; // NEW: path, job_idx, dir_idx

    // Separate vectors for sign files (no order number) with job indices
    std::vector<std::pair<std::string, int>> sign_files;
    std::vector<std::pair<std::string, int>> sign_block_std_files;
    std::vector<std::pair<std::string, int>> sign_seed_files;
    std::vector<std::pair<std::string, int>> sign_sigma_files;
    std::vector<std::pair<std::string, int>> sign_accept_ratio_files;
    std::vector<std::pair<std::string, int>> sign_direction_accept_ratio_files;
    std::vector<std::pair<std::string, int>> sign_tracked_arg_files;
    std::vector<std::tuple<std::string, int, int>> sign_histogram_files; // NEW: path, job_idx, dir_idx

    // Regular expressions to match file patterns and extract order numbers and job indices
    std::regex sign_pattern(R"(output_sign_(\d+)\.txt)");
    std::regex sign_block_pattern(R"(output_sign_block_std_(\d+)\.txt)");
    std::regex sign_seed_pattern(R"(output_sign_seed_(\d+)\.txt)");
    std::regex sign_sigma_pattern(R"(output_sign_sigma_(\d+)\.txt)");
    std::regex sign_accept_pattern(R"(output_sign_accept_ratio_(\d+)\.txt)");
    std::regex sign_dir_accept_pattern(R"(output_sign_direction_accept_ratios_(\d+)\.txt)");
    std::regex sign_tracked_arg_pattern(R"(output_sign_tracked_arg_(\d+)\.txt)");
    std::regex sign_hist_pattern(R"(output_sign_histogram_dir_(\d+)_(\d+)\.txt)"); // NEW

    std::regex grid_pattern(R"(output_grid_mode_\d+_order_(\d+)_(\d+)\.txt)");
    std::regex grid_block_pattern(R"(output_grid_mode_\d+_order_(\d+)_block_std_(\d+)\.txt)");
    std::regex grid_seed_pattern(R"(output_grid_mode_\d+_order_(\d+)_seed_(\d+)\.txt)");
    std::regex grid_sigma_pattern(R"(output_grid_mode_\d+_order_(\d+)_sigma_(\d+)\.txt)");
    std::regex grid_accept_pattern(R"(output_grid_mode_\d+_order_(\d+)_accept_ratio_(\d+)\.txt)");
    std::regex grid_dir_accept_pattern(R"(output_grid_mode_\d+_order_(\d+)_direction_accept_ratios_(\d+)\.txt)");
    std::regex grid_tracked_arg_pattern(R"(output_grid_mode_\d+_order_(\d+)_tracked_arg_(\d+)\.txt)");
    std::regex grid_hist_pattern(R"(output_grid_mode_\d+_order_(\d+)_histogram_dir_(\d+)_(\d+)\.txt)"); // NEW

    std::regex beta_pattern(R"(output_beta_mode_\d+_order_(\d+)_(\d+)\.txt)");
    std::regex beta_block_pattern(R"(output_beta_mode_\d+_order_(\d+)_block_std_(\d+)\.txt)");
    std::regex beta_seed_pattern(R"(output_beta_mode_\d+_order_(\d+)_seed_(\d+)\.txt)");
    std::regex beta_sigma_pattern(R"(output_beta_mode_\d+_order_(\d+)_sigma_(\d+)\.txt)");
    std::regex beta_accept_pattern(R"(output_beta_mode_\d+_order_(\d+)_accept_ratio_(\d+)\.txt)");
    std::regex beta_dir_accept_pattern(R"(output_beta_mode_\d+_order_(\d+)_direction_accept_ratios_(\d+)\.txt)");
    std::regex beta_tracked_arg_pattern(R"(output_beta_mode_\d+_order_(\d+)_tracked_arg_(\d+)\.txt)");
    std::regex beta_hist_pattern(R"(output_beta_mode_\d+_order_(\d+)_histogram_dir_(\d+)_(\d+)\.txt)"); // NEW

    std::regex mc_pattern(R"(output_MC_order_(\d+)_(\d+)\.txt)");
    std::regex mc_seed_pattern(R"(output_MC_order_(\d+)_seed_(\d+)\.txt)");

    // Set to track processed files
    std::set<std::string> files_to_delete;
    // Iterate through all files in save_path
    for (const auto& entry : std::filesystem::directory_iterator(save_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            std::string filename = entry.path().filename().string();
            std::string full_path = entry.path().string();

            std::smatch match;
            int order_num = -1;
            int job_index = -1;
            bool is_block_std = false;
            bool is_seed_file = false;
            bool is_sigma_file = false;
            bool is_accept_file = false;
            bool is_dir_accept_file = false;
            bool is_tracked_arg_file = false;
            bool is_histogram_file = false; // NEW
            bool is_sign_file = false;

            // Check each pattern and extract job index
            if (std::regex_match(filename, match, sign_pattern)) {
                job_index = std::stoi(match[1]);
                sign_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_block_pattern)) {
                job_index = std::stoi(match[1]);
                sign_block_std_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_seed_pattern)) {
                job_index = std::stoi(match[1]);
                sign_seed_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_sigma_pattern)) {
                job_index = std::stoi(match[1]);
                sign_sigma_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_accept_pattern)) {
                job_index = std::stoi(match[1]);
                sign_accept_ratio_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_dir_accept_pattern)) {
                job_index = std::stoi(match[1]);
                sign_direction_accept_ratio_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_tracked_arg_pattern)) {
                job_index = std::stoi(match[1]);
                sign_tracked_arg_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_hist_pattern)) { // NEW
                int dir_index = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                sign_histogram_files.emplace_back(full_path, job_index, dir_index);
                is_sign_file = true; // Still a sign file, but special kind
            } else if (std::regex_match(filename, match, grid_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
            } else if (std::regex_match(filename, match, grid_block_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_block_std = true;
            } else if (std::regex_match(filename, match, grid_seed_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_seed_file = true;
            } else if (std::regex_match(filename, match, grid_sigma_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_sigma_file = true;
            } else if (std::regex_match(filename, match, grid_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_accept_file = true;
            } else if (std::regex_match(filename, match, grid_dir_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_dir_accept_file = true;
            } else if (std::regex_match(filename, match, grid_tracked_arg_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_tracked_arg_file = true;
            } else if (std::regex_match(filename, match, grid_hist_pattern)) { // NEW
                order_num = std::stoi(match[1]);
                int dir_index = std::stoi(match[2]);
                job_index = std::stoi(match[3]);
                is_histogram_file = true;
            } else if (std::regex_match(filename, match, beta_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
            } else if (std::regex_match(filename, match, beta_block_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_block_std = true;
            } else if (std::regex_match(filename, match, beta_seed_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_seed_file = true;
            } else if (std::regex_match(filename, match, beta_sigma_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_sigma_file = true;
            } else if (std::regex_match(filename, match, beta_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_accept_file = true;
            } else if (std::regex_match(filename, match, beta_dir_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_dir_accept_file = true;
            } else if (std::regex_match(filename, match, beta_tracked_arg_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_tracked_arg_file = true;
            } else if (std::regex_match(filename, match, beta_hist_pattern)) { // NEW
                order_num = std::stoi(match[1]);
                int dir_index = std::stoi(match[2]);
                job_index = std::stoi(match[3]);
                is_histogram_file = true;
            } else if (std::regex_match(filename, match, mc_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
            } else if (std::regex_match(filename, match, mc_seed_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_seed_file = true;
            }
            // Add to appropriate map (only if not a sign file)
            if (!is_sign_file && order_num != -1 && job_index != -1) {
                if (is_histogram_file) {
                    int dir_index = std::stoi(match[2]); // Re-extract for safety
                    histogram_files_by_order[order_num].emplace_back(full_path, job_index, dir_index);
                } else if (is_tracked_arg_file) {
                    tracked_arg_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_seed_file) {
                    seed_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_block_std) {
                    block_std_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_sigma_file) {
                    sigma_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_accept_file) {
                    accept_ratio_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_dir_accept_file) {
                    direction_accept_ratio_files_by_order[order_num].push_back({full_path, job_index});
                } else {
                    regular_files_by_order[order_num].push_back({full_path, job_index});
                }
            }
        }
    }
    // Lambda to sort files by job index
    auto sort_by_job_index = [](auto& files) {
        std::sort(files.begin(), files.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });
    };

    // Lambda to sort histogram files by job index, then by direction index
    auto sort_hist_files = [](auto& files) {
        std::sort(files.begin(), files.end(),
                  [](const auto& a, const auto& b) {
                      if (std::get<1>(a) != std::get<1>(b)) {
                          return std::get<1>(a) < std::get<1>(b); // sort by job_index
                      }
                      return std::get<2>(a) < std::get<2>(b); // then by dir_index
                  });
    };
    // Combine sign files (without order number) with job index
    if (!sign_files.empty()) {
        sort_by_job_index(sign_files);
        std::string output_filename = std::string(save_path) + "/output_I.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign block std files (without order number) with job index
    if (!sign_block_std_files.empty()) {
        sort_by_job_index(sign_block_std_files);
        std::string output_filename = std::string(save_path) + "/output_I_block_std.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_block_std_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign seed files (without order number) with job index
    if (!sign_seed_files.empty()) {
        sort_by_job_index(sign_seed_files);
        std::string output_filename = std::string(save_path) + "/output_seeds.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_seed_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign sigma files (without order number) with job index
    if (!sign_sigma_files.empty()) {
        sort_by_job_index(sign_sigma_files);
        std::string output_filename = std::string(save_path) + "/output_sigmas.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_sigma_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign accept ratio files (without order number) with job index
    if (!sign_accept_ratio_files.empty()) {
        sort_by_job_index(sign_accept_ratio_files);
        std::string output_filename = std::string(save_path) + "/output_accept_ratios.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_accept_ratio_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign direction accept ratio files (without order number) with job index
    if (!sign_direction_accept_ratio_files.empty()) {
        sort_by_job_index(sign_direction_accept_ratio_files);
        std::string output_filename = std::string(save_path) + "/output_direction_accept_ratios.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_direction_accept_ratio_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign tracked arg files (without order number) with job index
    if (!sign_tracked_arg_files.empty()) {
        sort_by_job_index(sign_tracked_arg_files);
        std::string output_filename = std::string(save_path) + "/output_tracked_arg.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_tracked_arg_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign histogram files - NEW
    if (!sign_histogram_files.empty()) {
        sort_hist_files(sign_histogram_files);
        std::string output_filename = std::string(save_path) + "/output_histograms.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            output_file << "job_idx,direction_idx,bin_idx,count,probability" << std::endl;
            for (const auto& [input_file, job_idx, dir_idx] : sign_histogram_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    std::getline(input, line); // Skip header 1
                    std::getline(input, line); // Skip header 2
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << dir_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine regular files by order with job index
    for (auto& [order_num, files] : regular_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_I_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine block std files by order with job index
    for (auto& [order_num, files] : block_std_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_block_std_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine seed files by order with job index
    for (auto& [order_num, files] : seed_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_seeds_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sigma files by order with job index
    for (auto& [order_num, files] : sigma_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_sigmas_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine accept ratio files by order with job index
    for (auto& [order_num, files] : accept_ratio_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_accept_ratios_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine direction accept ratio files by order with job index
    for (auto& [order_num, files] : direction_accept_ratio_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_direction_accept_ratios_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine tracked arg files by order with job index
    for (auto& [order_num, files] : tracked_arg_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_tracked_arg_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine histogram files by order - NEW
    for (auto& [order_num, files] : histogram_files_by_order) {
        sort_hist_files(files);
        std::string output_filename = std::string(save_path) + "/output_histograms_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            output_file << "job_idx,direction_idx,bin_idx,count,probability" << std::endl;
            for (const auto& [input_file, job_idx, dir_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    std::getline(input, line); // Skip header 1
                    std::getline(input, line); // Skip header 2
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << dir_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }


    // Delete the original files
    for (const auto& file : files_to_delete) {
        try {
            std::filesystem::remove(file);
        } catch (const std::exception& e) {
            std::cerr << "Error deleting " << file << ": " << e.what() << std::endl;
        }
    }

    std::cout << "File combination complete." << std::endl;

}

inline void combine_order_update_files(int wait_seconds = 10) {
    // Wait for all files to be written
    std::cout << "Waiting " << wait_seconds << " seconds for all files to be written..." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(wait_seconds));
    std::cout << "Starting file combination for order update MCMC results..." << std::endl;

    // Convert save_path to string
    std::string save_dir = std::string(save_path);

    // Get order range from parameters
    const int order_min = N_MIN_ORDER;
    const int order_max = N_MAX_ORDER;

    // Maps to store files by run index
    std::map<int, std::vector<std::string>> visits_files;
    std::map<int, std::vector<std::string>> signs_files;
    std::map<int, std::vector<std::string>> seed_files;
    std::map<int, std::vector<std::string>> accept_files;
    std::map<int, std::vector<std::string>> params_files;
    std::map<int, std::vector<std::string>> sigma_files;
    std::map<int, std::vector<std::string>> summary_files;
    std::map<int, std::vector<std::string>> prelim_files;

    // Regular expressions to match file patterns
    std::regex visits_pattern(R"(visits_run_(\d+)\.txt)");
    std::regex signs_pattern(R"(signs_run_(\d+)\.txt)");
    std::regex seed_pattern(R"(seed_run_(\d+)\.txt)");
    std::regex accept_pattern(R"(accept_ratios_run_(\d+)\.txt)");
    std::regex params_pattern(R"(params_run_(\d+)\.txt)");
    std::regex sigma_pattern(R"(sigma_run_(\d+)\.txt)");
    std::regex summary_pattern(R"(summary_run_(\d+)\.txt)");
    std::regex prelim_pattern(R"(prelim_results_run_(\d+)\.txt)");

    // Set to track files to delete after combining
    std::set<std::string> files_to_delete;

    // Scan directory for files
    for (const auto& entry : std::filesystem::directory_iterator(save_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            std::string filename = entry.path().filename().string();
            std::string full_path = entry.path().string();
            std::smatch match;

            if (std::regex_match(filename, match, visits_pattern)) {
                int run_idx = std::stoi(match[1]);
                visits_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            } else if (std::regex_match(filename, match, signs_pattern)) {
                int run_idx = std::stoi(match[1]);
                signs_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            } else if (std::regex_match(filename, match, seed_pattern)) {
                int run_idx = std::stoi(match[1]);
                seed_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            } else if (std::regex_match(filename, match, accept_pattern)) {
                int run_idx = std::stoi(match[1]);
                accept_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            } else if (std::regex_match(filename, match, params_pattern)) {
                int run_idx = std::stoi(match[1]);
                params_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            } else if (std::regex_match(filename, match, sigma_pattern)) {
                int run_idx = std::stoi(match[1]);
                sigma_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            } else if (std::regex_match(filename, match, summary_pattern)) {
                int run_idx = std::stoi(match[1]);
                summary_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            } else if (std::regex_match(filename, match, prelim_pattern)) {
                int run_idx = std::stoi(match[1]);
                prelim_files[run_idx].push_back(full_path);
                files_to_delete.insert(full_path);
            }
        }
    }

    // Lambda to combine files of a specific type
    auto combine_files = [&save_dir, order_min, order_max](
        const std::map<int, std::vector<std::string>>& files_map,
        const std::string& output_suffix,
        bool include_run_index = true,
        const std::string& header = "") {

        if (files_map.empty()) return;

        std::string output_filename = save_dir + "/output_orders_range_" +
                                     std::to_string(order_min) + "_to_" +
                                     std::to_string(order_max) + "_" + output_suffix + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            if (!header.empty()) {
                output_file << header << std::endl;
            }

            // Process files in order of run index
            for (const auto& [run_idx, file_list] : files_map) {
                for (const auto& input_path : file_list) {
                    std::ifstream input(input_path);
                    if (input.is_open()) {
                        std::string line;
                        while (std::getline(input, line)) {
                            // Skip comment lines for some file types
                            if ((line.empty() || line[0] == '#')) {
                                if (!header.empty()) continue; // Skip headers if we have our own
                            }

                            if (include_run_index) {
                                output_file << run_idx << "," << line << std::endl;
                            } else {
                                output_file << line << std::endl;
                            }
                        }
                        input.close();
                    }
                }
            }
            output_file.close();
            std::cout << "Created: " << output_filename << std::endl;
        }
    };

    // Combine visits files
    combine_files(visits_files, "visits", true, "run_index,order,visits");

    // Combine signs files
    combine_files(signs_files, "signs", true, "run_index,order,mean_sign");

    // Combine seed files
    combine_files(seed_files, "seeds", true, "run_index,seed");

    // Combine acceptance ratio files
    combine_files(accept_files, "accept_ratios", false, "");

    // Combine parameter files
    combine_files(params_files, "parameters", false, "");

    // Combine sigma files
    combine_files(sigma_files, "sigmas", false, "");

    // Combine summary files (comprehensive output)
    if (!summary_files.empty()) {
        std::string output_filename = save_dir + "/output_orders_range_" +
                                     std::to_string(order_min) + "_to_" +
                                     std::to_string(order_max) + "_comprehensive_summary.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            output_file << "# Comprehensive Summary of All Runs\n";
            output_file << "# Order Range: [" << order_min << ", " << order_max << "]\n";
            output_file << "# Number of runs: " << summary_files.size() << "\n";
            output_file << "# " << std::string(70, '=') << "\n\n";

            for (const auto& [run_idx, file_list] : summary_files) {
                output_file << "# RUN INDEX: " << run_idx << "\n";
                output_file << "# " << std::string(70, '-') << "\n";

                for (const auto& input_path : file_list) {
                    std::ifstream input(input_path);
                    if (input.is_open()) {
                        std::string line;
                        while (std::getline(input, line)) {
                            output_file << line << std::endl;
                        }
                        input.close();
                    }
                }
                output_file << "\n";
            }
            output_file.close();
            std::cout << "Created: " << output_filename << std::endl;
        }
    }

    // Combine preliminary tuning results
    if (!prelim_files.empty()) {
        std::string output_filename = save_dir + "/output_orders_range_" +
                                     std::to_string(order_min) + "_to_" +
                                     std::to_string(order_max) + "_preliminary_tuning.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            output_file << "# Preliminary Parameter Tuning Results\n";
            output_file << "# run_index P_ADD P_REM accept_add accept_rem difference selected\n";

            for (const auto& [run_idx, file_list] : prelim_files) {
                for (const auto& input_path : file_list) {
                    std::ifstream input(input_path);
                    if (input.is_open()) {
                        std::string line;
                        bool found_selected = false;
                        std::string selected_params;

                        while (std::getline(input, line)) {
                            if (line.find("# Selected:") != std::string::npos) {
                                selected_params = line;
                                found_selected = true;
                            } else if (!line.empty() && line[0] != '#') {
                                output_file << run_idx << " " << line;
                                if (found_selected && line == selected_params.substr(12)) {
                                    output_file << " SELECTED";
                                }
                                output_file << std::endl;
                            }
                        }
                        input.close();
                    }
                }
            }
            output_file.close();
            std::cout << "Created: " << output_filename << std::endl;
        }
    }

    // Create a master index file
    {
        std::string index_filename = save_dir + "/output_orders_range_" +
                                    std::to_string(order_min) + "_to_" +
                                    std::to_string(order_max) + "_INDEX.txt";
        std::ofstream index_file(index_filename);

        if (index_file.is_open()) {
            index_file << "# Master Index of Combined Files\n";
            index_file << "# Order Range: [" << order_min << ", " << order_max << "]\n";
            index_file << "# Generated: " << std::chrono::system_clock::now().time_since_epoch().count() << "\n";
            index_file << "# Number of runs processed: " << std::max({
                visits_files.size(), signs_files.size(), seed_files.size(),
                accept_files.size(), params_files.size(), sigma_files.size()
            }) << "\n\n";

            index_file << "Files created:\n";
            index_file << "- visits.txt: Visit counts per order for each run\n";
            index_file << "- signs.txt: Mean signs per order for each run\n";
            index_file << "- seeds.txt: RNG seeds used for each run\n";
            index_file << "- accept_ratios.txt: Acceptance ratios and counts\n";
            index_file << "- parameters.txt: P_ADD, P_REM values used\n";
            index_file << "- sigmas.txt: Initial and final sigma values\n";
            index_file << "- comprehensive_summary.txt: Complete details of all runs\n";
            index_file << "- preliminary_tuning.txt: Parameter tuning results\n";

            index_file.close();
            std::cout << "Created index file: " << index_filename << std::endl;
        }
    }

    // Delete original files
    std::cout << "Cleaning up " << files_to_delete.size() << " original files..." << std::endl;
    for (const auto& file : files_to_delete) {
        try {
            std::filesystem::remove(file);
        } catch (const std::exception& e) {
            std::cerr << "Error deleting " << file << ": " << e.what() << std::endl;
        }
    }

    std::cout << "File combination complete for order range [" << order_min << ", " << order_max << "]" << std::endl;
}

#endif //FILE_COMBINATION_H
