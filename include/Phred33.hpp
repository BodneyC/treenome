/********************************************************************
 * Filename: Phred33.hpp [C++ header code]
 *
 * Description: File containing the probabilities of correct base
 *				reads as per Phred33 quality scores (1 -
 *P_error)
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
namespace GTH {
extern int score_sys;

const double phred_qualities[2][47] = {
    {0,        0.205672, 0.369043, 0.498813, 0.601893, 0.683772, 0.748811,
     0.800474, 0.841511, 0.874107, 0.9,      0.920567, 0.936904, 0.949881,
     0.960189, 0.968377, 0.974881, 0.980047, 0.984151, 0.987411, 0.99,
     0.992057, 0.99369,  0.994988, 0.996019, 0.996838, 0.997488, 0.998005,
     0.998415, 0.998741, 0.999,    0.999206, 0.999369, 0.999499, 0.999602,
     0.999684, 0.999749, 0.9998,   0.999842, 0.999874, 0.9999,   0.999921},
    {-2.16228, -1.51189, -0.995262, -0.584893, -0.258925, 0,        0.205672,
     0.369043, 0.498813, 0.601893,  0.683772,  0.748811,  0.800474, 0.841511,
     0.874107, 0.9,      0.920567,  0.936904,  0.949881,  0.960189, 0.968377,
     0.974881, 0.980047, 0.984151,  0.987411,  0.99,      0.992057, 0.99369,
     0.994988, 0.996019, 0.996838,  0.997488,  0.998005,  0.998415, 0.998741,
     0.999,    0.999206, 0.999369,  0.999499,  0.999602,  0.999684, 0.999749,
     0.9998,   0.999842, 0.999874,  0.9999,    0.999921}};
} // namespace GTH
