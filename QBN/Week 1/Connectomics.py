import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def closestDistanceBetweenLines(a0, a1, b0, b1, clampAll=False, clampA0=False, clampA1=False, clampB0=False,
                                clampB1=False):
    ''' Given two lines defined by numpy.array pairs (a0,a1,b0,b1)
        Return the closest points on each segment and their distance
    '''

    # If clampAll=True, set all clamps to True
    if clampAll:
        clampA0 = True
        clampA1 = True
        clampB0 = True
        clampB1 = True

    # Calculate denomitator
    A = a1 - a0
    B = b1 - b0
    magA = np.linalg.norm(A)
    magB = np.linalg.norm(B)

    _A = A / magA
    _B = B / magB

    cross = np.cross(_A, _B)
    denom = np.linalg.norm(cross) ** 2

    # If lines are parallel (denom=0) test if lines overlap.
    # If they don't overlap then there is a closest point solution.
    # If they do overlap, there are infinite closest positions, but there is a closest distance
    if not denom:
        d0 = np.dot(_A, (b0 - a0))

        # Overlap only possible with clamping
        if clampA0 or clampA1 or clampB0 or clampB1:
            d1 = np.dot(_A, (b1 - a0))

            # Is segment B before A?
            if d0 <= 0 >= d1:
                if clampA0 and clampB1:
                    if np.absolute(d0) < np.absolute(d1):
                        return a0, b0, np.linalg.norm(a0 - b0)
                    return a0, b1, np.linalg.norm(a0 - b1)


            # Is segment B after A?
            elif d0 >= magA <= d1:
                if clampA1 and clampB0:
                    if np.absolute(d0) < np.absolute(d1):
                        return a1, b0, np.linalg.norm(a1 - b0)
                    return a1, b1, np.linalg.norm(a1 - b1)

        # Segments overlap, return distance between parallel segments
        return None, None, np.linalg.norm(((d0 * _A) + a0) - b0)

    # Lines criss-cross: Calculate the projected closest points
    t = (b0 - a0)
    detA = np.linalg.det([t, _B, cross])
    detB = np.linalg.det([t, _A, cross])

    t0 = detA / denom
    t1 = detB / denom

    pA = a0 + (_A * t0)  # Projected closest point on segment A
    pB = b0 + (_B * t1)  # Projected closest point on segment B

    # Clamp projections
    if clampA0 or clampA1 or clampB0 or clampB1:
        if clampA0 and t0 < 0:
            pA = a0
        elif clampA1 and t0 > magA:
            pA = a1

        if clampB0 and t1 < 0:
            pB = b0
        elif clampB1 and t1 > magB:
            pB = b1

        # Clamp projection A
        if (clampA0 and t0 < 0) or (clampA1 and t0 > magA):
            dot = np.dot(_B, (pA - b0))
            if clampB0 and dot < 0:
                dot = 0
            elif clampB1 and dot > magB:
                dot = magB
            pB = b0 + (_B * dot)

        # Clamp projection B
        if (clampB0 and t1 < 0) or (clampB1 and t1 > magB):
            dot = np.dot(_A, (pB - a0))
            if clampA0 and dot < 0:
                dot = 0
            elif clampA1 and dot > magA:
                dot = magA
            pA = a0 + (_A * dot)

    return pA, pB, np.linalg.norm(pA - pB)

def point_in_cube(x):
    return x * (np.random.rand(3) - .5)

our_data = {'threshold':[],'connection probability':[],'number of connections':[]}
paper = {'threshold':[],'connection probability':[],'number of connections':[]}

for i in range(10):
    # simulation parameters
    cube_side = 10e-6  # [m]  # 1 we look at 10 μm sides for this
    threshold_dist =  1e-6 # [m] # 1 μm
    neurite_density = 10.61e18  # [m^(-3)] # density is 10.61 μm3. 
    num_samples = int(neurite_density * cube_side**3)

    # empty array to save distance samples in
    dist_samples = np.empty(num_samples)

    for i in tqdm(range(num_samples)):
        # generate line segments
        P1 = point_in_cube(cube_side)
        P2 = point_in_cube(cube_side)
        Q1 = point_in_cube(cube_side)
        Q2 = point_in_cube(cube_side)

        # calculate the closest distance between 2 lines
        P, Q, d = closestDistanceBetweenLines(P1, P2, Q1, Q2)

        # save distance
        dist_samples[i] = d

    # plot histogram
    our_data['threshold'].append(threshold_dist*10**6)
    our_data['connection probability'].append(np.sum(np.where(dist_samples<threshold_dist,1,0))/num_samples*100)
    our_data['number of connections'].append(np.sum(np.where(dist_samples<threshold_dist,1,0))/cube_side**3/10**9)

    threshold_dist = np.sort(dist_samples)[int(7.2e8* cube_side**3 * 10**9)] # such that Number of connections ~ 7.2*10^8 
    paper['threshold'].append(threshold_dist*10**6)
    paper['connection probability'].append(np.sum(np.where(dist_samples<threshold_dist,1,0))/num_samples*100)
    paper['number of connections'].append(np.sum(np.where(dist_samples<threshold_dist,1,0))/cube_side**3/10**9)

our_data['threshold'] = np.mean(our_data['threshold'])
our_data['connection probability'] = np.mean(our_data['connection probability'])
our_data['number of connections'] = np.mean(our_data['number of connections'])
paper['threshold'] = np.mean(paper['threshold'])
paper['connection probability'] = np.mean(paper['connection probability'])
paper['number of connections'] = np.mean(paper['number of connections'])

print('\n')
print('Averaged over 10 trials')
print('\n')
print('Our data:')
print('Threshold:',round(our_data['threshold'],2),'μm')
print('Connection probability:',str(round(our_data['connection probability'],3))+'%')
print('Number of connections in cubic millimeter: {:.1e}'.format(our_data['number of connections']))
print('\n')
print('Paper data:')
print('Threshold:',round(paper['threshold'],2),'μm')
print('Connection probability:',str(round(paper['connection probability'],3))+'%')
print('Number of connections in cubic millimeter: {:.1e}'.format(paper['number of connections']))

fig,ax = plt.subplots()
ax.hist(dist_samples,bins=100)
ax.set_xlabel('distance between neurites (m)')
ax.set_ylabel('number of trials')
ax.set_title('One trial of {} samples'.format(num_samples))
plt.show()