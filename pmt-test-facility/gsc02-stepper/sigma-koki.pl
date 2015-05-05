#!/usr/bin/perl -w
use GSC02;

use Sys::Hostname;
use strict;
use Device::SerialPort; 
use Socket;
use Time::HiRes qw(usleep nanosleep);

sub port_setup;
sub gsc_write;
sub hexwrite;

###################################################
# MAIN
###################################################
my $debug=0;
my $host = hostname();
## my $serial_device="/dev/ttyS0"; ## for pc with serial line
my $serial_device="/dev/ttyUSB0"; ## for pc with usb->serial adapter

my $controller=new GSC02($serial_device,1,2);



#TEST#

##$controller->goToHVOrigin();
##$controller->goToVOrigin();
##$controller->moveX("-",1);
##$controller->moveY("-",1);

$controller->move("w","+",10000);
#$controller->jog("1","+");
##$controller->move("1","-",1000);
##$controller->move("1","+",1000);



print($controller->isready(),"\n");
exit;


