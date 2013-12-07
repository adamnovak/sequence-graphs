import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter

schema = avro.schema.parse(open("schemas/AlleleGroup.avsc").read())

writer = DataFileWriter(open("allelegroups.avro", "w"), DatumWriter(), schema)

# Make a minimal AlleleGroup
ag = {
    "id": "ag1",
    "fivePrime": "ag1-5'",
    "threePrime": "ag1-3'",
    "contig": "chr1",
    "start": 0,
    "end": 10
    # No sequence, so we're the same as the reference.
}
    
    
writer.append(ag)
writer.close()

reader = DataFileReader(open("allelegroups.avro", "r"), DatumReader())
for group in reader:
    print group
reader.close()
