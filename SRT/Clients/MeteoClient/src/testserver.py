import socket
import sys

from socket import socket, AF_INET, SOCK_DGRAM
data = 'r ta01'
port = 5000
hostname = '192.168.202.66'
udp = socket(AF_INET,SOCK_DGRAM)
udp.sendto(data, (hostname, port))
received = udp.recv(1024)

print 
print received
