import bioCanon


def test_generate_new_codes():
    checking = "z"
    switch = {}
    groups = {'A': ['1', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              'B': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.1', 0, 0, 0, 0, 0, 0, 0],
              'C': ['z', 'z.1', 'z.1.1', 'z.1.1.z', 'z.1.1.z.1', 0, 0, 0, 0, 0, 0, 0],
              'D': ['z', 'z.z', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              'E': ['z', 'z.1', 'z.1.z', 'z.1.z.1', 'z.1.1.1.z', 0, 0, 0, 0, 0, 0, 0],
              'F': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.z', 0, 0, 0, 0, 0, 0],
              'G': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.z', 0,
                    0, 0, 0, 0],
              'H': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                    'z.1.1.1.z.1.1.z', 0, 0, 0, 0],
              'I': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                    'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.z', 0, 0, 0],
              'J': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                    'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.z', 0, 0],
              'K': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.z', 0, 0, 0, 0, 0, 0, 0],
              'L': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.z', 0, 0, 0, 0, 0, 0, 0],
              'M': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.1', 0, 0, 0, 0, 0, 0, 0],
              'N': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.1', 0, 0, 0, 0, 0, 0, 0],
              'O': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                    'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                    'z.1.1.1.z.1.1.1.1.1.z', 0],
              'P': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                    'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                    'z.1.1.1.z.1.1.1.1.1.z', 0],
              'Q': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                    'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                    'z.1.1.1.z.1.1.1.1.1.1', 0],
              'R': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                    'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                    'z.1.1.1.z.1.1.1.1.1.1', 0],
              'S': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.z', 0, 0, 0, 0, 0, 0],
              'T': ['z', 'z.1', 'z.1.1', 'z.1.1.z', 'z.1.1.z.z', 0, 0, 0, 0, 0, 0, 0],
              'U': ['z', 'z.1', 'z.1.1', 'z.1.1.z', 'z.1.1.z.z', 0, 0, 0, 0, 0, 0, 0],
              'V': ['z', 'z.z', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
    used = {1: ['1', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z', 'z',
                'z', 'z', 'z', 'z', 'z'],
            2: ['z.1', 'z.1', 'z.z', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1',
                'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.1', 'z.z'],
            3: ['z.1.1', 'z.1.1', 'z.1.z', 'z.1.1', 'z.1.1', 'z.1.1', 'z.1.1', 'z.1.1', 'z.1.z',
                'z.1.z', 'z.1.z', 'z.1.z', 'z.1.1', 'z.1.1', 'z.1.1', 'z.1.1', 'z.1.1', 'z.1.1',
                'z.1.1'],
            4: ['z.1.1.1', 'z.1.1.z', 'z.1.z.1', 'z.1.1.1', 'z.1.1.1', 'z.1.1.1', 'z.1.1.1',
                'z.1.1.1', 'z.1.z.z', 'z.1.z.z', 'z.1.z.z', 'z.1.z.z', 'z.1.1.1', 'z.1.1.1',
                'z.1.1.1', 'z.1.1.1', 'z.1.1.1', 'z.1.1.z', 'z.1.1.z'],
            5: ['z.1.1.1.1', 'z.1.1.z.1', 'z.1.1.1.z', 'z.1.1.1.z', 'z.1.1.1.z', 'z.1.1.1.z',
                'z.1.1.1.z', 'z.1.1.1.z', 'z.1.z.z.z', 'z.1.z.z.z', 'z.1.z.z.1', 'z.1.z.z.1',
                'z.1.1.1.z', 'z.1.1.1.z', 'z.1.1.1.z', 'z.1.1.1.z', 'z.1.1.1.z', 'z.1.1.z.z',
                'z.1.1.z.z'],
            6: ['z.1.1.1.z.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1', 'z.1.1.1.z.1', 'z.1.1.1.z.1',
                'z.1.1.1.z.1', 'z.1.1.1.z.1', 'z.1.1.1.z.1', 'z.1.1.1.z.1', 'z.1.1.1.z.z'],
            7: ['z.1.1.1.z.1.z', 'z.1.1.1.z.1.1', 'z.1.1.1.z.1.1', 'z.1.1.1.z.1.1', 'z.1.1.1.z.1.1',
                'z.1.1.1.z.1.1', 'z.1.1.1.z.1.1', 'z.1.1.1.z.1.1'],
            8: ['z.1.1.1.z.1.1.z', 'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1',
                'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1'],
            9: ['z.1.1.1.z.1.1.1.z', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1',
                'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1'],
            10: ['z.1.1.1.z.1.1.1.1.z', 'z.1.1.1.z.1.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                 'z.1.1.1.z.1.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1'],
            11: ['z.1.1.1.z.1.1.1.1.1.z', 'z.1.1.1.z.1.1.1.1.1.z', 'z.1.1.1.z.1.1.1.1.1.1',
                 'z.1.1.1.z.1.1.1.1.1.1']}
    iterators = [0, 'B']
    bioCanon.generate_new_codes(checking, switch, groups, used, iterators)
    compare = {'A': ['1', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               'B': ['2', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.1', 0, 0, 0, 0, 0, 0, 0],
               'C': ['z', 'z.1', 'z.1.1', 'z.1.1.z', 'z.1.1.z.1', 0, 0, 0, 0, 0, 0, 0],
               'D': ['z', 'z.z', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               'E': ['z', 'z.1', 'z.1.z', 'z.1.z.1', 'z.1.1.1.z', 0, 0, 0, 0, 0, 0, 0],
               'F': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.z', 0, 0, 0, 0, 0, 0],
               'G': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.z', 0,
                     0, 0, 0, 0],
               'H': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                     'z.1.1.1.z.1.1.z', 0, 0, 0, 0],
               'I': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                     'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.z', 0, 0, 0],
               'J': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                     'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.z', 0, 0],
               'K': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.z', 0, 0, 0, 0, 0, 0, 0],
               'L': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.z', 0, 0, 0, 0, 0, 0, 0],
               'M': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.1', 0, 0, 0, 0, 0, 0, 0],
               'N': ['z', 'z.1', 'z.1.z', 'z.1.z.z', 'z.1.z.z.1', 0, 0, 0, 0, 0, 0, 0],
               'O': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                     'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                     'z.1.1.1.z.1.1.1.1.1.z', 0],
               'P': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                     'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                     'z.1.1.1.z.1.1.1.1.1.z', 0],
               'Q': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                     'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                     'z.1.1.1.z.1.1.1.1.1.1', 0],
               'R': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.1', 'z.1.1.1.z.1.1',
                     'z.1.1.1.z.1.1.1', 'z.1.1.1.z.1.1.1.1', 'z.1.1.1.z.1.1.1.1.1',
                     'z.1.1.1.z.1.1.1.1.1.1', 0],
               'S': ['z', 'z.1', 'z.1.1', 'z.1.1.1', 'z.1.1.1.z', 'z.1.1.1.z.z', 0, 0, 0, 0, 0, 0],
               'T': ['z', 'z.1', 'z.1.1', 'z.1.1.z', 'z.1.1.z.z', 0, 0, 0, 0, 0, 0, 0],
               'U': ['z', 'z.1', 'z.1.1', 'z.1.1.z', 'z.1.1.z.z', 0, 0, 0, 0, 0, 0, 0],
               'V': ['z', 'z.z', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
    assert groups == compare
