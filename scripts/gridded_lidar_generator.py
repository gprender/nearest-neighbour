# gridded_lidar_generator.py
# A modified version of the random lidar generator script
# which creates a lidar file with a semi-regular grid distribution.

import random

def main():
    resolution = 4  # number of times to subdivide each axis
    bucketsize = 8  # number of points to generate in each cell
    xrange = (0,2**resolution)
    yrange = (0,2**resolution)
    zrange = (50,100)

    with open("grid.txt", 'w') as outfile:
        outfile.write(f'% min x y z          {xrange[0]} {yrange[0]} {zrange[0]}\n')
        outfile.write(f'% max x y z          {xrange[1]} {yrange[1]} {zrange[1]}\n')
        for i in range(2**resolution):
            for j in range(2**resolution):
                for k in range(bucketsize):
                    coord = (
                        round(random.uniform(i+0.01, i+0.99), 2),
                        round(random.uniform(j+0.01, j+0.99), 2),
                        round(random.uniform(zrange[0], zrange[1]), 2)
                    )
                    coord_str = [str(c) for c in coord]
                    for s in coord_str:
                        outfile.write(f'{s} ')
                    outfile.write('\n')

if __name__ == '__main__':
    main()