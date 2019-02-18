import sys
import socket

HOST = 'localhost'    #server name goes in here(Your IP)
PORT = 8001             	#port goes here
socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
socket.connect((HOST,PORT))
file = sys.argv[1]+"_testing.dat"			#name of tranformed sequence from SVMTriP
#$workdir.$filename . "_" . $epilength . "aa_prediction.txt"
outname = sys.argv[1]+"_"+sys.argv[2]+"aa_prediction.txt"	#name of prediction
with open( file, 'rb') as file_to_send:
    for data in file_to_send:
        socket.sendall(data)							#send the transformed sequence to the server
#print 'end'
#socket.close()

data = socket.recv(8001)				#gather the data of the prediction from the same socket
out = open(outname,"w")					#write to the output file expected from SVMTriP
out.write(data)

socket.close()
