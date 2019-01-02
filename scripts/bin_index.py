import sys

class interval(object):
    def __init__(self, start, end, *data):
        self.start = start
        self.end = end
        if not self.valid():
            raise RuntimeError("[interval] Error: interval start larger than"
                " end, [{},{}]".format(self.start,self.end))
        if len(data) > 0:
            self.data = data[0]
        else:
            self.data = None


    def valid(self):
        if self.start <= self.end:
            return 1
        else:
            return 0
    
    @property
    def size(self):
        return self.end - self.start + 1
    
    def overlap(self, other):
        if self.start > other.end or self.end < other.start:
            return 0
        return 1
    
    def __str__(self):
        return "{}\t{}".format(self.start, self.end)

class BinIndex(object):
    def __init__(self):
        self.db = dict()
        self.num_bins = 37449
        self.num_bin_levels = 6
        self._binFirstShift = 14
        self._binNextShift  = 3
        self._binOffsetsExtended = [0 for i in range(self.num_bin_levels)]
        # binOffsets = {4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};
        for i in range(self.num_bin_levels-2, -1, -1):
            self._binOffsetsExtended[i] = self._binOffsetsExtended[i+1] + (
                1 << ((self.num_bin_levels-2-i) * 3))
    
    def _getBin(self, start, end):
        end -= 1
        start >>= self._binFirstShift
        end >>= self._binFirstShift
        for i in range(self.num_bin_levels):
            if (start == end):
                return self._binOffsetsExtended[i] + start
            start >>= self._binNextShift
            end >>= self._binNextShift
        return -1
    
    def add_interval(self, _interval):
        bin_num = self._getBin(_interval.start, _interval.end)
        if (bin_num < 0 or bin_num > self.num_bins):
            print("[BinIndex] Error: Received illegal bin "
                "number {} from _getBin call.".format(bin_num))
        if bin_num not in self.db:
            self.db[bin_num] = [_interval]
        else:
            self.db[bin_num].append(_interval)
    
    def build_db(self, _intervals):
        for i in _intervals:
            self.add_interval(i)
    
    def get_overlap(self, _interval):
        result = []
        start_bin = _interval.start >> self._binFirstShift
        end_bin = _interval.end >> self._binFirstShift
        for i in range(self.num_bin_levels):
            offset = self._binOffsetsExtended[i]
            for j in range(start_bin+offset, end_bin+offset+1):
                if j not in self.db:
                    continue
                
                for k in self.db[j]:
                    if _interval.overlap(k):
                        result.append(k)
            start_bin >>= self._binNextShift
            end_bin >>= self._binNextShift
        if (len(result)>0):
            return result
        return 0

class bed_reader(object):
    def __init__(self, bed_fn):
        self.bed = bed_fn
    
    def __iter__(self):
        with open(self.bed, "r") as io:
            for line in io:
                fields = line.strip().split("\t")
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if start > end:
                    raise RuntimeError("Bad bed format start > end.{}".format(
                        line.strip()))
                if len(fields) > 3:
                    data = fields[3:]
                else:
                    data = 0
                if data:
                    _interval = interval(start, end, data)
                else:
                    _interval = interval(start, end)
                yield chrom, _interval

def build_db_from_bed(bed):
    bed_db = dict()
    bed_iter = bed_reader(bed)
    for chrom, _interval in bed_iter:
            if  chrom not in bed_db:
                bed_db[chrom] = BinIndex()
                bed_db[chrom].add_interval(_interval)
            else:
                bed_db[chrom].add_interval(_interval)
    return bed_db

def main():
    query_bed_iter = bed_reader(sys.argv[1])
    subject = build_db_from_bed(sys.argv[2])
    
    for chrom, _query_interval in query_bed_iter:
        if chrom in subject:
            _get = subject[chrom].get_overlap(_query_interval)
            if _get:
                for i in _get:
                    overlap = min(_query_interval.end, i.end) - max(_query_interval.start, i.start) + 1
                    if overlap/i.size >= 0.5 and overlap/_query_interval.size >= 0.5:
                        print(chrom+"\t",end="")
                        print(_query_interval,end="")
                        print("\t",end="")
                        print(chrom+"\t",end="")
                        print(i)

if __name__ == "__main__":
    main()
