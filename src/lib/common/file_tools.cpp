
#include <fstream>

namespace FileTools {

	bool is_file_exist(std::string fileName)
	{
		std::ifstream infile(fileName);
		return infile.good();
	}
}