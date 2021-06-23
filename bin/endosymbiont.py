import sys

def seq_length(fasta):
    total_length = 0
    for i in fasta:
        if(i[0] != '>'):
            total_length += len(i) - 1
    return total_length


def determine_endosymbiont(length, file):

    with open(length, 'r') as endosymbiont_seq_length:
        target_length = int(endosymbiont_seq_length.readline())
    
    with open(file) as contigs:
        contigs = contigs.readlines()
    
    query = []
    for i in range(0, len(contigs)):
        query.append(contigs[i])
        if seq_length(query) > 0.95 * target_length and seq_length(query) < 1.05 * target_length:
            with open('endosymbiont.fa', 'w') as endosymbiont:
                for item in query:
                    endosymbiont.write("%s" % item)
            break
    with open('endosymbiont.fa', 'w') as endosymbiont:
        for item in query:
            endosymbiont.write("%s" % item)
    pass

def main():
    length = sys.argv[1]
    file = sys.argv[2]
    determine_endosymbiont(length, file)
    pass

if __name__ == "__main__":
    main()