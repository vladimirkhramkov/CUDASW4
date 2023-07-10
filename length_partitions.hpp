#ifndef LENGTH_PARTITIONS_HPP
#define LENGTH_PARTITIONS_HPP

#include <array>
#include <limits>

//length k is in partition i if boundaries[i-1] < k <= boundaries[i]

constexpr auto getLengthPartitionBoundaries(){
    //const int l0 = 64, l1 = 128, l2 = 192, l3 = 256, l4 = 320, l5 = 384, l6 = 512, l7 = 640, l8 = 768, l9 = 1024, l10 = 5000, l11 = 50000;
	//const int l0 = 64, l1 = 128, l2 = 192, l3 = 256, l4 = 320, l5 = 384, l6 = 448, l7 = 512, l8 = 640, l9 = 768, l10 = 1024, l11 = 5000, l12 = 50000;
	//const int l0 = 64, l1 = 128, l2 = 192, l3 = 256, l4 = 320, l5 = 384, l6 = 448, l7 = 512, l8 = 640, l9 = 768, l10 = 896, l11 = 1024, l12 = 7000, l13 = 50000;
	
    //constexpr int numLengthPartitions = 15;
	//const int l0 = 64, l1 = 128, l2 = 192, l3 = 256, l4 = 320, l5 = 384, l6 = 448, l7 = 512, l8 = 576, l9 = 640, l10 = 768, l11 = 896, l12 = 1024, l13 = 7000, l14 = 50000;
	//std::array<int, numLengthPartitions> boundaries{l0,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14};

	constexpr int numLengthPartitions = 18;
	std::array<int, numLengthPartitions> boundaries{
		64,
		128,
		192,
		256,
		320,
		384,
		448,
		512,
		576,
		640,
		704,
		768,
		832,
		896,
		960,
		1024,
		7000,
		std::numeric_limits<int>::max()-1
	};

    return boundaries;
}
    



#endif