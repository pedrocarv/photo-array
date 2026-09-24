#include "Adafruit_Si7021.h"

Adafruit_Si7021 sensor = Adafruit_Si7021();
int command;

void setup() {
  Serial.begin(9600);

  // wait for serial port to open
  while (!Serial) {
    delay(10);
  }

  //Serial.println("Si7021 test!");
  
  if (!sensor.begin()) {
    //Serial.println("Did not find Si7021 sensor!");
    while (true)
      ;
  }
//
//  Serial.print("Found model ");
//  switch(sensor.getModel()) {
//    case SI_Engineering_Samples:
//      Serial.print("SI engineering samples"); break;
//    case SI_7013:
//      Serial.print("Si7013"); break;
//    case SI_7020:
//      Serial.print("Si7020"); break;
//    case SI_7021:
//      Serial.print("Si7021"); break;
//    case SI_UNKNOWN:
//    default:
//      Serial.print("Unknown");
//  }
  Serial.print(" Rev(");
  Serial.print(sensor.getRevision());
  Serial.print(")");
  Serial.print(" Serial #"); Serial.print(sensor.sernum_a, HEX); Serial.println(sensor.sernum_b, HEX);
}

char recvOneChar() {
 if (Serial.available() > 0) {
 return Serial.read();
 }
}

void loop() {
 if (Serial.available() > 0) {
  command = Serial.read();
  if (command == 't') {
  Serial.println(sensor.readTemperature(), 2);
 }
 if (command == 'h'){
   Serial.println(sensor.readHumidity(), 2);
 }
}
}
//
//void loop() {
//  //Serial.print("Humidity:    ");
//  Serial.print(sensor.readHumidity(), 2);
//  //Serial.print("\tTemperature: ");
//  Serial.print(sensor.readTemperature(), 2);
//  delay(1000);
//}
