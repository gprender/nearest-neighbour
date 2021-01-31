# random_lidar_generator.py
# A quick little script that generates a lidar file in .txt format,
# with output similar to that of LAStools' las2txt. 

import random

def main():
	num_points = 100000
	xrange = (0.01, 499.99)
	yrange = (0.01, 499.99)
	zrange = (50,100)
	
	with open("rand.txt", 'w') as outfile:
		outfile.write(f'% min x y z          {xrange[0]} {yrange[0]} {zrange[0]}\n')
		outfile.write(f'% max x y z          {xrange[1]} {yrange[1]} {zrange[1]}\n')
		for i in range(num_points):
			coord = (
				round(random.uniform(xrange[0], xrange[1]), 2),
				round(random.uniform(yrange[0], yrange[1]), 2),
				round(random.uniform(zrange[0], zrange[1]), 2)
			)
			coord_str = [str(c) for c in coord]
			for s in coord_str:
				outfile.write(f'{s} ')
			outfile.write('\n')

if __name__ == '__main__':
    main()