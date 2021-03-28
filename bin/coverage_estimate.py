alignment_rate = open('alignment_rate.txt', 'r')               # Get alignment rate from bowtie2 output
line = alignment_rate.readlines()[-1]
rate = ''
for i in line:
    if i == '%': break
    rate += i
rate = float(rate)/100
alignment_rate.close()

host_count = open('host_count.txt', 'r')            # Get no. of bases in reference genome
count = int(host_count.readline())
host_count.close()

base_count = open('base_count.txt', 'r')            # Get no. of bases in read files
bases = int(base_count.readline()) * 2

coverage = (rate * bases) / count                   # Calculate coverage

coverage_file = open('coverage.txt', 'w')
coverage_file.write(str(coverage))
coverage_file.close
