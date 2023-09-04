

def compute_barrier_props(barriers_mod: dict, barriers_names: dict, barrier_name: str) -> dict:
    """ compute barrier geometry
    """
    barrier_props = {}

    # height/depth
    barrier_h_d = get_barrier_height_and_depth(barriers_mod, barriers_names, barrier_name)
    barrier_props.update(barrier_h_d)
    
    # radius
    barrier_r = get_barrier_radius(barriers_mod, barrier_name)
    barrier_props.update(barrier_r)

    return barrier_props

def get_barrier_height_and_depth(barriers_mod: dict, barriers_names: dict, barrier_name: str) -> dict:
    '''Get height (m) top (msl) nd bottom (msl) of the barrier   barrier_name.
        As it is now only the main barrier is used as input, but the code here
        is general so one might e.g. loop over all barriers if needed
    '''
    top    = []
    bottom = []
    for bname in barriers_names[barrier_name]:                       #Loop over the sections of this barrier
        top.append(barriers_mod['top_msl'][bname])
        bottom.append(barriers_mod['bottom_msl'][bname])

    barrier_props = {}

    barrier_props['height'] = max(bottom) - min(top)
    barrier_props['top']    = min(top)
    barrier_props['bottom'] = max(bottom)

    return barrier_props

def get_barrier_radius(barriers_mod: dict, barriers_names: dict, barrier_name: str):
    '''Height-averaged radius if the barrier varies in diameter'''

    heights = []
    diams   = []
    
    #Collect diameters and heights
    for bname in barriers_names[barrier_name]:
        heights.append(barriers_mod['bottom_msl'][bname] - barriers_mod['top_msl'][bname])
        diams.append(barriers_mod['diameter_m'][bname])
        
    #Do the avaraging
    avg_diam = 0
    for diam, height in zip(diams, heights):
        avg_diam += diam*height
    avg_diam /= sum(heights)

    barrier_props = {}
    
    barrier_props['radius'] = avg_diam/2.0

    return barrier_props
