def get_distinct_clones( handle, with_count=False ):

    ## LOOKS THROUGH TEXT FILE OF CLASSIFIERS AND WRITES NEW FILE CONTAINING ALL DISTINCT CLASSIFIERS, OPTIONALLY WITH COUNT OF ALL DISTINCT CLASSIFIERS
    ## with_count=True writes file with counts of all distinct classifiers
    
    from string import Template
    import collections as coll
    from operator import itemgetter, attrgetter

    write_to = open("distinct_clones.txt", "w")

    if with_count:
        stemplate = Template('$count $element')
        d = coll.defaultdict(int)
        for line in handle:
            elements = line.rstrip("\n")
            d[elements] += 1
        d_sorted = sorted(d.items(), key=itemgetter(1), reverse=True)
        for k in d_sorted:
            f_seq = stemplate.substitute( count = k[1], element = k[0])
            print >> write_to, f_seq
    else:
        stemplate = Template('$element')
        d = coll.defaultdict(int)
        for line in handle:
            elements = line.rstrip("\n")
            if elements not in d:
                d[elements] = 1
                f_seq = stemplate.substitute(element = elements)
                print >> write_to, f_seq
    
    handle.close()
    write_to.close()