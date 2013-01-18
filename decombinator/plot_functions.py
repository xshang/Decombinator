def plot_v_usage( handle, savefilename="Vusage", order="frequency"):

    ## PLOTS V GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    tags = open("tags_trbv.txt", "rU")
    num_genes = 0
    for line in tags:
        num_genes += 1

    freq_vector_v = [0]*num_genes
    for line in handle:
        elements = line.rstrip("\n")
        freq_vector_v[int(string.split(elements)[0])] += 1

    if order=="frequency":        
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_v)
        percent_usage_v = [0]*num_genes
        for i in range(num_genes):
            percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/V12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27-1', 'V28-1', 'V29-1', 'V3-1', 'V30-1', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9')
        v_linked = [0]*len(percent_usage_v)
        for i in range(len(percent_usage_v)):
            v_linked[i] = (gene_list_v[i], percent_usage_v[i])
        sorted_v = sorted(v_linked, key=itemgetter(1))
        v_labels = [0]*len(sorted_v)
        v_percents = [0]*len(sorted_v)
        for j in range(len(sorted_v)):
            v_labels[j] = sorted_v[j][0]
            v_percents[j] = sorted_v[j][1]
        pos_v = np.arange(num_genes)+ 1
        plt.figure()
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.yticks( pos_v, v_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)
        
    elif order=="chromosome":
        total = sum(freq_vector_v)
        fv = [0]*num_genes
        for i in range(num_genes):
            fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        gene_list_v = ('10-1', '10-2', '10-3', '11-1', '11-2', '11-3', '12-3/4', '12-5', '13', '14', '15', '16', '18', '19', '2', '20-1', '24-1', '25-1', '27-1', '28-1', '29-1', '3-1', '30-1', '4-1', '4-2', '4-3', '5-1', '5-4', '5-5', '5-6', '5-8', '6-1', '6-4', '6-5', '6-6', '6-8', '6-9', '7-2', '7-3', '7-4', '7-6', '7-7', '7-8', '7-9', '9')
        chromosome_order = [ 14, 21, 23, 26, 31, 24, 25, 37, 32, 38, 44, 0, 3, 1, 4, 33, 39, 27, 34, 28, 40, 29, 35, 41, 36, 42, 30, 43, 8, 2, 5, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 22 ]
        gene_list_v = [ gene_list_v[i] for i in chromosome_order]
        fv = [fv[i] for i in chromosome_order]

        ind = np.arange(num_genes)
        width = 0.25

        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects1 = ax.bar(ind, fv, width, color='yellow')

        ax.set_ylabel('Frequency', fontsize = 10)
        ax.set_xticks(ind+1*width)
        ax.set_xticklabels(gene_list_v)
        plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 10)
        plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 10)
        plt.grid(True)

        plt.savefig(str(savefilename)+'.png', dpi=300)
        
    handle.close()
    tags.close()

def plot_j_usage( handle, savefilename="Jusage", order="frequency"):

    ## PLOTS J GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    tags = open("tags_trbj.txt", "rU")
    num_genes = 0
    for line in tags:
        num_genes += 1

    freq_vector_j = [0]*num_genes
    for line in handle:
        elements = line.rstrip("\n")
        freq_vector_j[int(string.split(elements)[1])] += 1

    if order=="frequency":    
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_j)
        percent_usage_j = [0]*num_genes
        for i in range(num_genes):
            percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
        j_linked = [0]*len(percent_usage_j)
        for i in range(len(percent_usage_j)):
            j_linked[i] = (gene_list_j[i], percent_usage_j[i])
        sorted_j = sorted(j_linked, key=itemgetter(1))
        j_labels = [0]*len(sorted_j)
        j_percents = [0]*len(sorted_j)
        for j in range(len(sorted_j)):
            j_labels[j] = sorted_j[j][0]
            j_percents[j] = sorted_j[j][1]
        pos_j = np.arange(num_genes)+ 1
        plt.figure()
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.yticks( pos_j, j_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)

    elif order=="chromosome":
        total = sum(freq_vector_j)
        fj = [0]*num_genes
        for i in range(num_genes):
            fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        gene_list_j = ('1-1', '1-2', '1-3', '1-4', '1-5', '1-6', '2-1', '2-2', '2-3', '2-4', '2-5', '2-6', '2-7')

        ind = np.arange(num_genes)
        width = 0.25

        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects1 = ax.bar(ind, fj, width, color='red')

        ax.set_ylabel('Frequency', fontsize = 10)
        ax.set_xticks(ind+width)
        ax.set_xticklabels(gene_list_j)
        plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 10)
        plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 10)
        plt.grid(True)

        plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags.close()

def plot_del_v( handle, savefilename="Vdels"):

    ## PLOTS V GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_v = [0]*50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_v[int(string.split(elements)[2])] += 1

    total = sum(deletions_v)
    for i in range(len(deletions_v)):
        deletions_v[i] = deletions_v[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_v, width, color='yellow')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of V germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()

def plot_del_j( handle, savefilename="Jdels"):

    ## PLOTS J GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_j = [0]*50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_j[int(string.split(elements)[3])] += 1

    total = sum(deletions_j)
    for i in range(len(deletions_j)):
        deletions_j[i] = deletions_j[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_j, width, color='red')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of J germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()

def plot_vj_joint_dist( handle, savefilename="VJusage", tags_v = open("tags_trbv.txt", "rU"), tags_j = open("tags_trbj.txt", "rU")):

    ## PLOTS VJ JOINT GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec

    num_v = 0
    for line in tags_v:
        num_v += 1

    num_j = 0
    for line in tags_j:
        num_j += 1
        
    joint_distribution = np.zeros((num_v,num_j))
    for line in handle:
        elements = line.rstrip("\n")

        v = int(string.split(elements)[0])
        j = int(string.split(elements)[1])

        joint_distribution[v,j] += 1

    joint_distribution = joint_distribution / sum(sum(joint_distribution))
    gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/V12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27-1', 'V28-1', 'V29-1', 'V3-1', 'V30-1', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9')
    gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')

    pos_v = np.arange(num_v)+ 1
    pos_j = np.arange(num_j)+ 1
    
    plt.figure()
    plt.pcolor(joint_distribution)
    pos_ticks_v = pos_v-0.5
    pos_ticks_j = pos_j-0.5
    plt.yticks( pos_ticks_v, gene_list_v)
    plt.xticks( pos_ticks_j, gene_list_j)
    plt.colorbar()
    plt.pcolor(joint_distribution)
    yticklabels = plt.getp(plt.gca(), 'yticklabels')
    plt.setp(yticklabels, fontsize='8')
    xticklabels = plt.getp(plt.gca(), 'xticklabels')
    plt.setp(xticklabels, fontsize='8')
    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags_v.close()
    tags_j.close()

def plot_insert_lengths( handle, savefilename="InsertLengths" ):

    ## PLOTS DISTRIBUTION OF INSERT LENGTH BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    maxi = 100
    insert_lengths = [0]*maxi

    handle = open("male1_total_unique.txt", "rU")
    for line in handle:
        elements = line.rstrip("\n")

        classifier = string.split(elements)
        if len(classifier) == 5:
            insert_lengths[len(classifier[4])] += 1
        else:
            insert_lengths[0] += 1

    total = sum(insert_lengths)
    for i in range(len(insert_lengths)):
        insert_lengths[i] = insert_lengths[i] / float(total)

    ind = np.arange(maxi)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, insert_lengths, width, color='blue')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of nucleotides', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.xlim((0,50))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
