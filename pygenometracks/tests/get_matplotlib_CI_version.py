def get_CI_mpl_version():
    with open('requirements_CI.txt', 'r') as f:
        for line in f:
            if line.startswith('matplotlib =='):
                # Remove potential comments
                if '#' in line:
                    line = line.split('#')[0]
                # Remove potential spaces
                line = line.strip()
                return(line.replace('matplotlib ==', ''))
