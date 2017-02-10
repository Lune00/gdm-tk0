#include "io.hpp"

using namespace std;

int main(int argc,char **argv)
{
	Sample * spl = new Sample;
	sample_read(argv[1], *spl);
	sample_write_tapioK("packing.cin", *spl);

	delete spl;
	return 0;
}
