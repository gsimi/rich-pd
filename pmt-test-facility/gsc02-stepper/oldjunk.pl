




#open and setup the serial port
my $port = port_setup($serial_device);

#write something
printf ("doing something\n");
$port->lookclear; 

gsc_write($port,"Q:\r\n");
##-      990,         0,X,K,R
my ($count,$line)=$port->read(27); # will read _up to_ 255 chars
printf("$line\n");

gsc_write($port,"M:1+P10000\r\n");
gsc_write($port,"G\r\n");



sub port_setup($){
    my $serial_device=shift;
    
    my $port = Device::SerialPort->new($serial_device);
    defined ($port) || die "can't open port $serial_device";
    $port->baudrate(9600);
    $port->parity("none");  ## none, odd, even
    $port->handshake("rts");## none, rts, xoff
    $port->databits(8);     ## 5,6,7,8
    $port->stopbits(1);     ## 1,2
    $port->read_char_time(0);
    $port->read_const_time(1);
    $port->write_settings ;
    return $port;
}

sub hexwrite($$){
    my $port=shift;
    my $hc=shift;
    $port->write(chr(hex("$hc\r\n")));
}

sub gsc_write($$){
    my $port=shift;
    my $output_string=shift;
    $output_string=$output_string."\r\n"; ## Add CR+LF
    
    my $count_out = $port->write($output_string);
    warn "write failed\n"         unless ($count_out);
    warn "write incomplete, written $count_out, expected".length($output_string)."\n" unless ($count_out eq length($output_string));
    usleep(100000); ##wait 0.1 sec

}
