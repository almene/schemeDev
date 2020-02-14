from bioCanon import bioCanon


def test_id_conflict():
    required_tiles = [359, 515, 1895, 9797, 9864]
    flanking = 15
    all_variable= [68, 254, 305, 359, 366, 507, 515, 605, 684, 816, 909, 945, 954, 957, 963, 966,
                   1017, 1092, 1142, 1180, 1200, 1325, 1366, 1486, 1628, 1673, 1679, 1697, 1893,
                   1895, 2025, 2049, 2050, 2139, 2310, 2379, 2527, 2566, 2610, 2685, 2866, 2981,
                   3115, 3349, 3538, 3547, 4801, 5021, 5063, 5065, 5066, 5067, 5068, 5110, 5371,
                   5392, 5806, 5889, 5930, 5931, 5932, 5976, 5979, 6032, 6040, 6085, 6136, 6190,
                   6519, 6523, 6529, 6625, 6651, 6699, 6736, 6777, 6837, 6991, 7207, 7258, 7266,
                   7288, 7480, 7554, 7915, 8075, 8086, 8197, 8209, 8368, 8398, 8479, 8856, 9121,
                   9174, 9222, 9312, 9323, 9327, 9417, 9633, 9797, 9864, 10664]
    case = bioCanon.id_conflict(required_tiles, flanking, all_variable)
    compare = {359: [366], 515: [507], 1895: [1893]}
    assert case == compare
