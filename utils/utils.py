
def gen_dict_extract(key, var):
    if hasattr(var, 'iteritems'):
        for k, v in var.iteritems():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result

def find_dict_key(key, dictionary):
    for k, v in dictionary.items():
        if k == key:
            return v
        if isinstance(v, dict):
            res = find_dict_key(key, v)
            if res != 0:
                return res
        if isinstance(v, list):
            for d in v:
                res = find_dict_key(key, d)
                if res != 0:
                    return res

    return 0
