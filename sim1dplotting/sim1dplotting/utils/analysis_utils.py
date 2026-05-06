def map_dx(Nx, dx):
    
    if dx >= 1:
        return dx
    
    mapping = {
        627: 800,
        802: 625,
        202: 2500,
        1002: 500,
        2002: 250,
        4002: 125,
        5002: 100,
        8002: 62.5,
    }
    
    return mapping.get(Nx, dx)
    