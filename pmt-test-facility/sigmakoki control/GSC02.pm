#!/usr/bin/perl
package GSC02;

use strict;
use Device::SerialPort;
use Socket;
use Time::HiRes qw(usleep nanosleep);
use Switch;
use warnings;

# class method
sub new() {
	## class name as first arg
	my $this = shift;
	my $class = ref($this) || $this;
	## allocates a hash and return a reference to it
	## data members are initialized in initialize()
	my $self = {
		_serial_device =>,
		_hstage        =>,
		_vstage        =>,
		_port          =>,
		_usleep        => 100000    ## 0.1sec;
	};

	## tells to the $self reference to has that is a GSC02
	bless $self, $class;
	$self->initialize(@_);
	return $self;
}

# instalnce methods
sub initialize() {
	my $self = shift;
	$self->{_serial_device} = shift;
	printf("Serial device $self->{_serial_device}\n");
	$self->{_hstage} = shift;
	printf( "Horizontal stage on axis %d\n", $self->{_hstage} );
	$self->{_vstage} = shift;
	printf( "Vertical stage on axis %d\n", $self->{_vstage} );

	my $port = Device::SerialPort->new( $self->{_serial_device} );
	defined($port) || die "can't open port $self->_serial_device";
	$self->{_port} = $port;
	$port->baudrate(9600);
	$port->parity("none");      ## none, odd, even
	$port->handshake("rts");    ## none, rts, xoff
	$port->databits(8);         ## 5,6,7,8
	$port->stopbits(1);         ## 1,2
	$port->read_char_time(0);
	$port->read_const_time(1);
	$port->write_settings;
	printf("Port configured\n");
}

#write a string on serial line
sub write() {
	my ( $self, $output_string ) = @_;
	$output_string = $output_string . "\r\n";    ## Add CR+LF
	my $count_out = $self->{_port}->write($output_string);
	warn "write failed\n" unless ($count_out);
	warn "write incomplete, written $count_out, expected"
	  . length($output_string) . "\n"
	  unless ( $count_out eq length($output_string) );
	usleep( $self->{_usleep} );                  ##wait 0.1 sec
}

sub checkaxis () {
	my ( $self, $axis ) = @_;
	unless (
		defined($axis)
		&& (   ( $axis eq "W" )
			|| ( $axis eq "w" )
			|| ( $axis == $self->{_hstage} )
			|| ( $axis == $self->{_vstage} ) )
	  )
	{
		do {
			print(
"Error : wrong axis. You have to specify one of those axes. W, 1 or 2\n"
			);
			chomp( $axis = <> );
		  } while ( ( $axis ne "W" )
			&& ( $axis ne "w" )
			&& ( $axis ne $self->{_hstage} )
			&& ( $axis ne $self->{_vstage} ) );
		return $axis;
	}
	return $axis;
}

sub checkdir() {
	my ( $self, $dir ) = @_;
	unless (
		defined($dir)
		&& (   ( $dir eq "+" )
			|| ( $dir eq "-" )
			|| ( $dir eq "++" )
			|| ( $dir eq "--" )
			|| ( $dir eq "-+" )
			|| ( $dir eq "+-" ) )
	  )
	{
		do {
			print(
"Error : wrong direction. You have to specify a direction. +,-,++,--,-+ or +-\n"
			);
			chomp( $dir = <> );
		  } while ( ( $dir ne "+" )
			&& ( $dir ne "-" )
			&& ( $dir ne "++" )
			&& ( $dir ne "--" )
			&& ( $dir ne "-+" )
			&& ( $dir ne "+-" ) );
		return $dir;
	}
	return $dir;
}

########################
#      goToOrigin      #
########################

sub goToOrigin() {
	my ( $self, $axis, $dir ) = @_;
	$self->waitready();
	$axis = $self->checkaxis($axis);
	$dir  = $self->checkdir($dir);
	if ( ( $axis eq "W" ) || ( $axis eq "w" ) ) {
		printf("Sending both axes to $dir origin\n");
		if ( ( $dir eq "+" ) || ( $dir eq "-" ) ) {
			$self->write("H:W${dir}${$dir}");
		}
		else {
			$self->write("H:W${dir}");
		}
		$self->write("G");
	}
	else {
		if (   $dir eq "++"
			|| $dir eq "+-"
			|| $dir eq "--"
			|| $dir eq "-+" )
		{
			do {
				printf(
"Error, direction $dir is not valid in this case. Enter + or -\n"
				);
				chomp( $dir = <> );
			  } while ( ( $dir ne ("-") )
				&& ( $dir ne ("+") ) );
		}
		if ( $axis == 1 ) {
			printf("Sending axis $axis to $dir origin\n");
			$self->write("H:${axis}${dir}");
			$self->write("G");
		}
		else {
			printf("Sending axis $axis to $dir origin\n");
			$self->write("H:${axis}${dir}");
			$self->write("G");
		}
	}
}

sub goToHOrigin() {
	my ( $self, $dir ) = @_;
	my $axis = $self->{_hstage};
	$self->goToOrigin( $axis, $dir );
}

sub goToVOrigin() {
	my ( $self, $dir ) = @_;
	my $axis = $self->{_vstage};
	$self->goToOrigin( $axis, $dir );
}

sub goToBothOrigin() {
	my ( $self, $dir1, $dir2 ) = @_;

	#	my $axis1 = $self->{_hstage};
	#	my $axis2 = $self->{_vstage};( $dir ne "+" )
	$self->goToOrigin( "W", $dir1, $dir2 );
}

########################
#         move         #
########################

sub move() {
	my ( $self, $axis, $dir, $steps1, $steps2 ) = @_;
	$self->waitready();
	$axis = $self->checkaxis($axis);
	$dir  = $self->checkdir($dir);
	if ( ( $axis eq ("W") ) || ( $axis eq ("w") ) ) {
		if ( !defined $steps1 ) {
			do {
				print(
					"Error : you have to specify a number of steps for axis 1\n"
				);
				chomp( $steps1 = <> );
			} while ( !defined $steps1 );
		}
		if ( ( $steps1 < 0 ) || ( $steps1 > 16777214 ) ) {
			do {
				printf(
"Error: number of steps $steps1 out of range. Enter a value between 0 and 16777214\n"
				);
				chomp( $steps1 = <> );
			} while ( ( $steps1 < 0 ) || ( $steps1 > 16777214 ) );
		}
		if ( !defined $steps2 ) {
			do {
				print(
					"Error : you have to specify a number of steps for axis 2\n"
				);
				chomp( $steps2 = <> );
			} while ( !defined $steps2 );
		}
		if ( ( $steps2 < 0 ) || ( $steps2 > 16777214 ) ) {
			do {
				printf(
"Error: number of steps $steps2 out of range. Enter a value between 0 and 16777214\n"
				);
				chomp( $steps2 = <> );
			} while ( ( $steps2 < 0 ) || ( $steps2 > 16777214 ) );
		}
		switch ($dir) {
			case "+" {
				$self->write("M:W+P${steps1}+P${steps2}");
				$self->write("G");
			}
			case "-" {
				$self->write("M:W-P${steps1}-P${steps2}");
				$self->write("G");
			}
			case "++" {
				$self->write("M:W+P${steps1}+P${steps2}");
				$self->write("G");
			}
			case "+-" {
				$self->write("M:W+P${steps1}-P${steps2}");
				$self->write("G");
			}
			case "-+" {
				$self->write("M:W-P${steps1}+P${steps2}");
				$self->write("G");
			}
			case "--" {
				$self->write("M:W-P${steps1}-P${steps2}");
				$self->write("G");
			}
		}
	}
	else {
		if ( !defined $steps1 ) {
			do {
				print(
					"Error : you have to specify a number of steps for axis \n"
				);
				chomp( $steps1 = <> );
			} while ( !defined $steps1 );
		}
		if ( ( $steps1 < 0 ) || ( $steps1 > 16777214 ) ) {
			do {
				printf(
"Error: number of steps $steps1 out of range. Enter a value between 0 and 16777214\n"
				);
				chomp( $steps1 = <> );
			} while ( ( $steps1 < 0 ) || ( $steps1 > 16777214 ) );
		}
		switch ($dir) {
			case "+" {
				$self->write("M:${axis}+P${steps1}");
				$self->write("G");
			}
			case "-" {
				$self->write("M:${axis}-P${steps1}");
				$self->write("G");
			}
			case "++" {
				$self->write("M:${axis}+P${steps1}");
				$self->write("G");
			}
			case "-+" {
				do {
					printf(
"Error, direction $dir is not valid in this case. Enter + or -\n"
					);
					chomp( $dir = <> );
				  } while ( ( $dir ne ("-") )
					&& ( $dir ne ("+") ) );
				$self->write("M:${axis}${dir}P${steps1}");
				$self->write("G");
			}
			case "+-" {
				do {
					printf(
"Error, direction $dir is not valid in this case. Enter + or -\n"
					);
					chomp( $dir = <> );
				  } while ( ( $dir ne ("-") )
					&& ( $dir ne ("+") ) );
				$self->write("M:${axis}${dir}P${steps1}");
				$self->write("G");
			}
			case "--" {
				$self->write("M:${axis}-P${steps1}");
				$self->write("G");
			}
		}
	}
}

sub moveBoth() {
	my ( $self, $dir, $steps1, $steps2 ) = @_;
	$self->move( "W", $dir, $steps1, $steps2 );
}

sub moveH() {
	my ( $self, $dir, $steps ) = @_;
	my $axis = $self->{_hstage};
	$self->move( $axis, $dir, $steps );
}

sub moveV() {
	my ( $self, $dir, $steps ) = @_;
	my $axis = $self->{_vstage};
	$self->move( $axis, $dir, $steps );
}

sub isint() {
	my ( $self, $value ) = @_;
	if ( $value =~ /^-?\d+\z/ ) {
		return 1;
	}
	return 0;
}

sub moveHMMeter() {
	my ( $self, $dir, $mm ) = @_;
	if ( !defined $mm || $mm < 0 || $mm > 100) {
		do {
			print("Error : you have to specify a number of mm for horizontal axis between 0 and 100 \n");
			chomp( $mm = <> );
		} while ( !defined $mm || ($mm < 0 || $mm > 100) );
	}
	my $steps = $mm*500;
	unless($self->isint($steps) && $steps>=0 && $steps<=50000){
		do {
			print("Error : ${mm}mm not correspond to a good number of steps. \n");
			chomp( $mm = <> );
			$steps=$mm*500;
		} while ( ! $self->isint($steps) ||  $steps<0 || $steps>50000);
	}
	$self->moveH( $dir, $steps );
}

sub moveVMMeter() {
	my ( $self, $dir, $mm ) = @_;
	if ( !defined $mm || $mm < 0 || $mm > 100) {
		do {
			print("Error : you have to specify a number of mm for vertical axis between 0 and 100 \n");
			chomp( $mm = <> );
		} while ( !defined $mm || ($mm < 0 || $mm > 100) );
	}
	my $steps = $mm*500;
	unless($self->isint($steps) && $steps>=0 && $steps<=50000){
		do {
			print("Error : ${mm}mm not correspond to a good number of steps. \n");
			chomp( $mm = <> );
			$steps=$mm*500;
		} while ( ! $self->isint($steps) ||  $steps<0 || $steps>50000);
	}
	$self->moveV( $dir, $steps );
}

########################
#   setLogicalOrigin   #
########################

sub setLogicalOrigin () {
	my ( $self, $axis ) = @_;
	if ( $self->checkaxis($axis) ) {
		$self->write("R:${axis}");
	}
}

sub setHLogicalOrigin () {
	my ($self) = @_;
	my $axis = $self->{_hstage};
	$self->setLogicalOrigin($axis);
	if ( $self->getVPosition() == 0 ) {
		print("The logical origin of vertical axis is set \n");
	}
	else {
		print("Error !");
		exit;
	}
}

sub setVLogicalOrigin () {
	my ($self) = @_;
	my $axis = $self->{_vstage};
	$self->setLogicalOrigin($axis);

	#	my $debug = $self->getVPosition();
	#	print("Debug : VPostion : $debug\n");
	if ( $self->getVPosition() == 0 ) {
		print("The logical origin of vertical axis is set \n");
	}
	else {
		print("Error !");
		exit;
	}
}

sub setBothLogicalOrigin() {
	my ($self) = @_;
	$self->setLogicalOrigin("W");
	if ( $self->getHPosition() == 0 && $self->getVPosition() == 0 ) {
		print("The logical origin of both axes is set \n");
	}
	else {
		print("Error !");
		exit;
	}
}

sub setSpeed() {
	my ( $self, $range, $minSpd1, $maxSpd1, $accelerationTime1, $minSpd2,
		$maxSpd2, $accelerationTime2 )
	  = @_;
	$self->waitready();
	unless ( $range == 1 || $range == 2 ) {
		printf("Error, non existing range $range\n");
		exit;
	}
	unless ( ( $range == 1 ) && ( $minSpd1 >= 1 ) && ( $minSpd1 <= 200 ) ) {
		printf("Error : speed $minSpd1 out of range\n");
		exit;
	}
	unless ( ( $range == 1 ) && ( $maxSpd1 >= 1 ) && ( $maxSpd1 <= 200 ) ) {
		printf("Error : speed $maxSpd1 out of range\n");
		exit;
	}
	unless ( ( $range == 2 )
		&& ( $minSpd2 >= 50 )
		&& ( $minSpd2 <= 20000 ) )
	{
		printf("Error : speed $minSpd2 out of range\n");
		exit;
	}
	unless ( ( $range == 2 )
		&& ( $maxSpd2 >= 50 )
		&& ( $maxSpd2 <= 20000 ) )
	{
		printf("Error : speed $maxSpd2 out of range\n");
		exit;
	}
	unless ( $accelerationTime1 <= 0 || $accelerationTime1 >= 1000 ) {
		printf("Error : acceleration time  $accelerationTime1 out of range\n");
		exit;
	}
	unless ( $accelerationTime2 <= 0 || $accelerationTime2 >= 1000 ) {
		printf("Error : acceleration time  $accelerationTime2 out of range\n");
		exit;
	}
	$self->write(
"D:${range}S${minSpd1}F${maxSpd1}R${accelerationTime1}S${minSpd2}F${maxSpd2}R${accelerationTime2}"
	);
}

sub freeMotor { }

sub holdMotor { }

########################
#     getPosition      #
########################

sub getPositions() {
	my ($self) = @_;
	$self->write("Q:");
	my ( $count, $buffer ) = $self->{_port}->read(127);
	chomp($buffer);

	#	print ("debug: position string is $buffer\n");
	#	print ("debug: position count is $count\n");
	my ( $axis1Position, $axis2Position ) = split( ',', $buffer );

	#	print("debug axis 1 pos :$axis1Position/n");
	$axis1Position =~ s/\s+//;
	$axis2Position =~ s/\s+//;
	my $hposition;
	my $vposition;
	if ( $self->{_hstage} == 1 ) {
		$hposition = $axis1Position;
		$vposition = $axis2Position;
	}
	else {
		$hposition = $axis2Position;
		$vposition = $axis1Position;
	}
	my @coordinates = ( $hposition, $vposition );

	#	print("my coordinates are : ($hposition,$vposition) \n");
	return @coordinates;
}

sub getHPosition() {
	my ($self) = @_;
	my $hPosition = ( $self->getPositions() )[0];
	print("Horizontal Position : $hPosition \n");
	return $hPosition;
}

sub getVPosition() {
	my ($self) = @_;
	my $vPosition = ( $self->getPositions() )[1];
	print("Vertical Position : $vPosition \n");
	return $vPosition;
}

sub getBothPosition() {
	my ($self) = @_;
	my @coordinates = $self->getPositions();
	print("Coordinates are : ($coordinates[0],$coordinates[1]) \n");
	return @coordinates;
}

sub jog() {
	my ( $self, $axis, $dir ) = @_;
	$self->waitready();
	$self->checkaxis($axis);
	$self->checkdir($dir);
	$self->write("J:${axis}${dir}");
	$self->stop( $self, $axis );
	$self->write("G");
}

sub stop() {
	my ( $self, $axis ) = @_;
	$self->checkaxis($axis);
	$self->write("L:${axis}");
}

sub emergencyStop () {
	my ($self) = @_;
	$self->write("L:E");
}

sub isready($) {
	my ($self) = @_;
	$self->write("!:");
	my ( $count, $buffer ) = $self->{_port}->read(3);
	chomp($buffer);
	if ( $buffer =~ "^R" ) {
		return 1;
	}
	return 0;
}

sub waitready($) {
	my ($self) = @_;
	while ( !$self->isready() ) { usleep(500000); }
}

sub version { }

1;
