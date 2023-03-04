def copy_dicts(dicts):
    new_dicts = []
    for dictionary in dicts:
        new_dicts.append(dictionary.copy())
    return new_dicts


def modify_dicts(dicts, modifier_list):
    new_dicts = copy_dicts(dicts)
    for modifier in modifier_list:
        print(modifier)
