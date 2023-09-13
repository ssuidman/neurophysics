import os

DB = 'db79'
# Update file paths based on the location of your MATLAB code
MATLAB_DIR = '/Users/sam/Documents/MATLAB/Code/db79'  # Path to your MATLAB code directory

pumcu_db_version_path = os.path.join(MATLAB_DIR, 'pumcu_db_version.txt')
pumcu_diagsIdName_path = os.path.join(MATLAB_DIR, 'pumcu_diagsIdName.txt')
pumcu_testsIdName_path = os.path.join(MATLAB_DIR, 'pumcu_testsIdName.txt')
pumcu_externalsIdName_path = os.path.join(MATLAB_DIR, 'pumcu_externalsIdName.txt')
pumcu_diagprevs_path = os.path.join(MATLAB_DIR, 'pumcu_diagprevs.txt')
pumcu_agemult_path = os.path.join(MATLAB_DIR, 'pumcu_agemult.txt')
pumcu_gendermult_path = os.path.join(MATLAB_DIR, 'pumcu_gendermult.txt')
MULT = 'improved'
maxage = 100

# Set up data structures
diagn = {}
test = {}
prev = {}
agemult = {}
malemult = {}
femalemult = {}

# Helper function to remove overlapping age mult intervals
def remove_overlapping_intervals(low, high, mult):
    non_overlapping_low = []
    non_overlapping_high = []
    non_overlapping_mult = []

    for i in range(len(low)):
        overlap = False
        for j in range(len(low)):
            if i != j:
                if high[i] <= low[j] or high[j] <= low[i]:
                    overlap = False
                else:
                    overlap = True
                    break

        if not overlap:
            non_overlapping_low.append(low[i])
            non_overlapping_high.append(high[i])
            non_overlapping_mult.append(mult[i])

    return non_overlapping_low, non_overlapping_high, non_overlapping_mult

# Read pumcu_db_version.txt
with open(pumcu_db_version_path, 'r') as fid:
    db_version = fid.readline().strip()
    print(f'database version: {db_version}')

# Read pumcu_diagsIdName.txt
with open(pumcu_diagsIdName_path, 'r') as fid:
    a = fid.readline().strip()
    while a:
        j = a.find('@')
        d = int(a[:j])
        diagn[d] = a
        a = fid.readline().strip()

Ndiags = len(diagn)
print(f'Ndiags = {Ndiags}')

# Read pumcu_testsIdName.txt
with open(pumcu_testsIdName_path, 'r') as fid:
    a = fid.readline().strip()
    while a:
        j1 = a.find('@')
        t = int(a[:j1])
        j2 = a.find('@', j1 + 1)
        test[t] = {'gender': a[j1 + 1:j2], 'name': a[j2 + 1:]}
        a = fid.readline().strip()

Ntests = len(test)
print(f'Ntests = {Ntests}')

# Read pumcu_externalsIdName.txt
with open(pumcu_externalsIdName_path, 'r') as fid:
    a = fid.readline().strip()
    while a:
        j1 = a.find('@')
        t = int(a[:j1])
        test[t]['name'] = a[j1 + 1:]
        test[t]['type'] = 6
        a = fid.readline().strip()

Ntests = len(test)
print(f'Ntests = {Ntests}')

# Read pumcu_diagprevs.txt
with open(pumcu_diagprevs_path, 'r') as fid:
    a = fid.readline()
    while a:
        j = a.find('@')
        d = int(a[:j])
        prev[d] = float(a[j + 1:])
        a = fid.readline()

# Read pumcu_agemult.txt
pumcu_agemult = []

with open(pumcu_agemult_path, 'r') as fid:
    for line in fid:
        pumcu_agemult.append(list(map(float, line.split())))

nagemult = [0] * Ndiags

for i in range(len(pumcu_agemult)):
    d = int(pumcu_agemult[i][0])
    nagemult[d] += 1
    if d not in agemult:
        agemult[d] = {'low': [], 'high': [], 'mult': []}

    agemult[d]['low'].append(pumcu_agemult[i][1])
    agemult[d]['high'].append(pumcu_agemult[i][2])
    agemult[d]['mult'].append(pumcu_agemult[i][3])

for d in range(Ndiags):
    if nagemult[d] == 0:
        agemult[d] = {'low': [], 'high': [], 'mult': []}

# Remove overlapping age mult intervals
maxmult = 0
minmult = 1000

for d in range(Ndiags):
    if nagemult[d] > 1:
        low = agemult[d]['low']
        high = agemult[d]['high']
        mult = agemult[d]['mult']
        maxmult = max(max(mult), maxmult)
        minmult = min(min(mult), minmult)
        non_overlapping_low, non_overlapping_high, non_overlapping_mult = remove_overlapping_intervals(low, high, mult)
        nagemult[d] = len(non_overlapping_low)
        agemult[d]['low'] = non_overlapping_low
        agemult[d]['high'] = non_overlapping_high
        agemult[d]['mult'] = non_overlapping_mult

# Switch statement equivalent using a dictionary
switch_dict = {
    'improved': {
        'description': 'adjust prevs to agemult',
        'normalrange': maxage,
        'prevcorr': None
    }
}

if MULT in switch_dict:
    print(switch_dict[MULT]['description'])
    print('ASSUME ALL AGES EQUALLY LIKELY')
    for d in range(Ndiags):
        if nagemult[d] > 0:
            low = agemult[d]['low']
            high = [min(x, maxage) for x in agemult[d]['high']]
            mult = agemult[d]['mult']
            normalrange = maxage - sum(high[i] - low[i] for i in range(len(low)))
            prevcorr = (normalrange + sum((high[i] - low[i]) * mult[i] for i in range(len(low)))) / maxage
            prev[d] = prev[d] / prevcorr

# Read pumcu_gendermult.txt
pumcu_gendermult = []

with open(pumcu_gendermult_path, 'r') as fid:
    for line in fid:
        pumcu_gendermult.append(list(map(int, line.split())))

for i in range(len(pumcu_gendermult)):
    d = int(pumcu_gendermult[i][0])
    m = pumcu_gendermult[i][2]

    if MULT == 'current':
        if pumcu_gendermult[i][1] == 0:
            femalemult[d] = m
        else:
            malemult[d] = m
    elif MULT == 'improved':
        if pumcu_gendermult[i][1] == 0:
            malemult[d] = 2 / (1 + m)
            femalemult[d] = 2 * m / (1 + m)
        else:
            malemult[d] = 2 * m / (1 + m)
            femalemult[d] = 2 / (1 + m)

# Adjust gender mults from diagsIdName
for d in range(Ndiags):
    a = diagn[d]
    j1 = a.find('@')
    j2 = a.find('@', j1 + 1)
    gender = a[j1 + 1:j2]
    if gender == 'female':
        malemult[d] = 0
    elif gender == 'male':
        femalemult[d] = 0
    diagn[d] = a[j2 + 1:]

print('GENDER INFORMATION FROM UMCU_TESTSIdNAME.TXT IGNORED')

# Read normal distributions (pumcu_normal.txt) - You can implement this part if needed
# Load pumcu_normal.txt and process it as required
