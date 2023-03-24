# ifndef IO_CPP
# define IO_CPP

# include <string>
# include <vector>

bool str2bool(std::string s);

void read_binary_float(std::string path, float * array, int n);

void write_binary_float(std::string path, float *array, int n);

void read_text_file(std::string path, std::vector<std::string> &elements);

std::string catch_parameter(std::string target, std::string file);

std::vector<std::string> split(std::string s, char delimiter);

# endif