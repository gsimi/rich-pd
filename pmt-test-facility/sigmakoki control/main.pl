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
my $debug = 0;
my $host  = hostname();
## my $serial_device="/dev/ttyS0"; ## for pc with serial line
my $serial_device = "/dev/ttyUSB0";    ## for pc with usb->serial adapter

my $controller = new GSC02( $serial_device, 1, 2 );

#TEST#
#
#printf($controller->isready());
#$controller->write("Q:");
#my ( $count, $buffer )=(0,"");
#print("count : $count /n");
#print("buffer : $buffer");
#$controller->write("G");
#print("issued move command\n");
#my ( $count, $buffer )=(0,"");
#while (!($buffer =~ "^R")) {
#		print ("writing !: on gsc\n");
#	$controller->write("!:");
#	( $count, $buffer ) =  $controller->{_port}->read(3);
#	print ("read buffer $buffer\n");
#print("read count $count\n");	
#		usleep(500000);
#  } 
#print("gsc is ready\n");
#$controller->write("Q:");
#( $count, $buffer ) = $controller->{_port}->read(100);
#print("count $count\n");
#print("buffer $buffer\n");


#$controller->setBothLogicalOrigin();

#print($controller->checkdir());


#$controller->goToOrigin();

#$controller->waitready();
#$controller->getHPosition();
#$controller->setBothLogicalOrigin();
#$controller->getHPosition();

$controller->moveBoth();
#$controller->waitready();
#print($controller->getPositions());
#$controller->setBothLogicalOrigin();
#$controller->getPositions();

#$controller->moveV("+",10000000);
#print(isint(0));

#print($controller->isint(01));

#$controller->moveVMMeter();
#$controller->getBothPosition();
#$controller->setBothLogicalOrigin();

#$controller->goToBothOrigin();
#$controller->setVLogicalOrigin();

#$controller->getPositions();
#$controller->getVPosition();

#$controller->goToVOrigin("+");
#$controller->getPosition();

#$controller->getPosition();

#$controller->checkdir("l");
#$controller->checkaxis("1");

#$controller->goToOrigin("W","+","+");
#$controller->goToOrigin("W","++");
#$controller->emergencyStop();

#$controller->moveH("+",1000000);
#$controller->moveV("-",10000);

#print($controller->isready(),"\n");
exit;

