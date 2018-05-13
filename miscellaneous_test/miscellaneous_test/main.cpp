#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>

using namespace::std;

void main()
{
	for (int i = 0; i < 10; i++)
	{
		ofstream file;
		file.open("test_writting.txt",fstream::app);
		file << to_string(i) << "\n";
		file.close();
	}
}