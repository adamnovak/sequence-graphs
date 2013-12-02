import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter

schema = avro.schema.parse(open("user.avsc").read())

writer = DataFileWriter(open("users.avro", "w"), DatumWriter(), schema)

alyssa = {"name": "Alyssa", "favorite_number": 256}
ben = {"name": "Ben", "favorite_number": 7, "favorite_color": "red"}

alyssa["favorite_user"] = ben
ben["favorite_user"] = alyssa

writer.append(alyssa)
writer.append(ben)
writer.close()

reader = DataFileReader(open("users.avro", "r"), DatumReader())
for user in reader:
    print user
reader.close()
