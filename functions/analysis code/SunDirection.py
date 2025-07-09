def sun_direction(UTdays, UTsecs, UTnsecs):
    import numpy as np
    #Esta funcion no recibe un array de numeros! solo valores numericos y no listas.
    # Constantes
    labtwist_degrees = -49.58
    pi = np.pi
    twopi = 2 * pi
    to_rad = pi / 180.0

    # Offset en días: diferencia entre el tiempo cero de SNO+ (1 Jan 2010) y Jan 0, 2000
    days_offset = 2451543.5 - 2455197.5

    # Días desde Jan 0, 2000
    days = UTdays + UTsecs / 86400.0 + UTnsecs / (1.0e9 * 86400.0) - days_offset

    # Constantes astronómicas
    ecl = (23.4393 - 3.563E-7 * days) * to_rad
    w = (282.9404 + 4.70935E-5 * days) * to_rad
    e = 0.016709 - 1.151E-9 * days
    ma = (356.0470 + 0.9856002585 * days) * to_rad

    # Eccentric anomaly
    ea = ma + e * np.sin(ma) * (1.0 + e * np.cos(ma))

    # Dirección hacia el Sol (coordenadas eclípticas)
    xv = np.cos(ea) - e
    yv = np.sqrt(1.0 - e * e) * np.sin(ea)

    # True anomaly
    v = np.arctan2(yv, xv)

    # Longitud solar
    sun_longitude = v + w

    # Coordenadas geocéntricas eclípticas
    xs = np.cos(sun_longitude)
    ys = np.sin(sun_longitude)
    vec1 = np.array([
        ys * np.cos(ecl),
        ys * np.sin(ecl),
        xs
    ])

    # Rotación terrestre diaria
    k0 = 0.27499
    k1 = 1.0 + 1.0 / 365.2425
    spin = k0 + k1 * days
    spin_angle = (spin - int(spin)) * twopi

    # Latitud y longitud del detector (en radianes)
    longitude = (81 + 12 / 60 + 5 / 3600) * to_rad
    latitude = (46 + 28 / 60 + 31.0 / 3600) * to_rad
    labtwist = labtwist_degrees * to_rad

    # Rotaciones: Z(labtwist) * X(latitude) * Y(longitude - spin_angle)
    def rotation_matrix_y(angle):
        return np.array([
            [np.cos(angle), 0, np.sin(angle)],
            [0, 1, 0],
            [-np.sin(angle), 0, np.cos(angle)]
        ])

    def rotation_matrix_x(angle):
        return np.array([
            [1, 0, 0],
            [0, np.cos(angle), -np.sin(angle)],
            [0, np.sin(angle), np.cos(angle)]
        ])

    def rotation_matrix_z(angle):
        return np.array([
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1]
        ])

    rot = (
        rotation_matrix_z(labtwist) @
        rotation_matrix_x(latitude) @
        rotation_matrix_y(longitude - spin_angle)
    )

    sun_dir = rot @ vec1  # Producto matriz * vector

    return sun_dir  # Es un array de 3 elementos: [x, y, z]