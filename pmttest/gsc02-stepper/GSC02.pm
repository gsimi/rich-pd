#!/usr/bin/perl
package GSC02;

use strict;
use Device::SerialPort; 
use Socket;
use Time::HiRes qw(usleep nanosleep);

# class method
sub new
{
    ## class name as first arg
    my $this = shift;
    my $class = ref($this) || $this;
    ## allocates a hash and return a reference to it
    ## data members are initialized in initialize()
    my $self = {    _serial_device => ,
		    _hstage        => ,
		    _vstage        => ,
		    _port          => ,
		    _usleep        => 100000 ## 0.1sec;
    };
    
    ## tells to the $self reference to has that is a GSC02
    bless $self, $class;
    $self->initialize(@_);
    return $self;
}

# instalnce methods
sub initialize{
    my $self=shift;
    $self->{_serial_device} = shift;
    printf ("Serial device $self->{_serial_device}\n");
    $self->{_hstage}        = shift;
    printf ("Horizontal stage on axis %d\n",$self->{_hstage});
    $self->{_vstage}        = shift;
    printf ("Vertical stage on axis %d\n",$self->{_vstage});

    my $port = Device::SerialPort->new($self->{_serial_device});
    defined ($port) || die "can't open port $self->_serial_device";
    $self->{_port} = $port;
    $port->baudrate(9600);
    $port->parity("none");  ## none, odd, even
    $port->handshake("rts");## none, rts, xoff
    $port->databits(8);     ## 5,6,7,8
    $port->stopbits(1);     ## 1,2
    $port->read_char_time(0);
    $port->read_const_time(1);
    $port->write_settings ;
    printf ("Port configured\n");
}

#write a string on serial line
sub write{
    my ($self,$output_string)=@_;
    $output_string=$output_string."\r\n"; ## Add CR+LF
    
    my $count_out = $self->{_port}->write($output_string);
    warn "write failed\n"         unless ($count_out);
    warn "write incomplete, written $count_out, expected".length($output_string)."\n" unless ($count_out eq length($output_string));
    usleep($self->{_usleep}); ##wait 0.1 sec
}

#sub checkaxis{
#    my ($self,$axis)=@_;
#    unless (($axis eq "W") or ($axis eq "w") or ($axis eq $self->{_hstage}) or ($axis eq $self->{_vstage}){
#	printf("Error, non existing axis $axis\n");
#	exit;
#   }
#}


sub checkaxis{
    my ($self,$axis)=@_;
   if (($axis eq "W") || ($axis eq "w")){
   }
	else{
  	  if (($axis != $self->{_hstage}) && ($axis != $self->{_vstage})) {
		printf("Error, non existing axis $axis\n");
		exit;
  	  }
	}
}

sub checkdir{
    my ($self,$dir)=@_;
    if (($dir =~/"+"/) && ($dir =~/"-"/)){
	printf("Error, non existing direction $dir\n");	
	exit;
    }
}

## goes to origin
sub goToOrigin($$$){
    my ($self, $axis, $dir) = @_;
    $self->waitready();
    print "Sending axis $axis to $dir origin\n";
    $self->write("H:${axis}${dir}");
    $self->write("G");
}
sub goToHOrigin($$) {
    my $self = shift;
    my $dir = shift || "+";
    checkdir($dir);
    my $axis=$self->{_hstage};
    $self->goToOrigin($axis,$dir);
}    

sub goToVOrigin($$) {
    my $self = shift;
    my $dir = shift || "+";
    checkdir($dir);
    my $axis=$self->{_vstage};
    $self->goToOrigin($axis,$dir);
}    

sub goToHVOrigin($$) {
    my $self = shift;
    my $dir = shift || "++";
    my $axis="W";
    $self->goToOrigin($axis,$dir);
}    

sub move($$$$){
    my ($self, $axis, $dir,$steps) = @_;
    $self->waitready();
    $self->checkaxis($axis);
    $self->checkdir($dir);
    printf("Moving axis $axis, direction $dir, steps $steps\n");
    if (($steps<0) || ($steps>16777214)){
	printf("Error: number of steps $steps out of range\n");
	exit;
    }
    if($axis eq ("W") or $axis eq ("w")){
	$self->write("M:W${dir}P${steps}${dir}P${steps}");
        $self->write("G");
    }
    else{
        $self->write("M:${axis}${dir}P${steps}");
        $self->write("G");
    }
}

#sub move($$$$$){
#    my ($self, $dir1,$steps1,$dir2,$steps2) = @_;
#    $self->waitready();
#    $self->checkaxis($axis);
#    $self->checkdir($dir);
#    printf("Moving axis $axis1, direction $dir1, steps $steps1\n");
#    if (($steps1<0) || ($steps1>16777214)){
#	printf("Error: number of steps $steps1 out of range\n");
#	exit;
#    }
#    if (($steps2<0) || ($steps2>16777214)){
#	printf("Error: number of steps $steps2 out of range\n");
#	exit;
#    }
#    if($axis eq ("W" || "w")){
#	$self->write("M:W${dir1}P${steps1}${dir2}P${steps2}");
#        $self->write("G");
#    }
#}



sub moveX($$$){    
    my ($self, $dir,$steps) = @_;    
    my $stage=$self->{_hstage};
    $self->move($stage,$dir,$steps);
}
sub moveY($$$){
    my ($self, $dir,$steps) = @_;    
    printf("direction $dir, steps $steps\n");
    $self->move($self->{_vstage},$dir,$steps);
}


sub jog($$$$){
 my ($self, $axis, $dir) = @_;
    $self->waitready();
    $self->checkaxis($axis);
    $self->checkdir($dir);
    $self->write("J:${axis}${dir}");
    $self->stop($self,$axis);
    $self->write("G");
}


sub stop($$$){
    my ($self, $axis) =@_;
    $self->checkaxis($axis);
    $self->write("L:${axis}");
}

sub emergencyStop{}
sub resetOrigin{}
sub setSpeed{}
sub freeMotor{}
sub holdMotor{}
sub getPosition{}

sub isready($){
    my ($self)=@_;
    $self->write("!:");
    my ($count,$buffer)=$self->{_port}->read(1);
    chomp($buffer);
    if ($buffer eq "R"){
	return 1;}
    return 0;
}

sub waitready($){
    my ($self)=@_;
    while (! $self->isready()) {usleep(500000);}
}

sub version{}

1;
