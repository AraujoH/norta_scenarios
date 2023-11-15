using AWSS3

# Set the credentials
aws_access_key_id = "YOUR_AWS_ACCESS_KEY_ID"
aws_secret_access_key = "YOUR_AWS_SECRET_ACCESS_KEY"

# Set the bucket and file path
bucket = "arpa-e-perform"
file_path = "path/to/your/file.hdf5"

# Download the file
s3_read(aws_access_key_id, aws_secret_access_key, bucket, file_path)


url = "https://data.openei.org/files/341/BA_load_day-ahead_fcst_2018.h5"
using HDF5
# Open the file and read its contents
f = h5open(url, "r")
data = read(f["/forecast"])
close(f)

# Show the contents of the file
show(data)


#https://data.openei.org/s3_viewer?bucket=arpa-e-perform&prefix=ERCOT%2F2018%2FLoad%2FActuals%2FBA_level%2F
