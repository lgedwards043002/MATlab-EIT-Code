#include <Wire.h>

// MUX select pins
const int S0 = 2;
const int S1 = 3;
const int S2 = 4;
const int S3 = 5;

struct Electrode {
  uint8_t muxID;
  uint8_t channel;
  uint8_t adsAddr;
};

// E2–E23
Electrode electrodes[22] = {
  {0, 0, 0x48}, {0, 1, 0x48}, {0, 2, 0x48}, {0, 3, 0x48},
  {0, 4, 0x48}, {0, 5, 0x48}, {0, 6, 0x48}, {0, 7, 0x48},
  {1, 0, 0x49}, {1, 1, 0x49}, {1, 2, 0x49}, {1, 3, 0x49},
  {1, 4, 0x49}, {1, 5, 0x49}, {1, 6, 0x49}, {1, 7, 0x49},
  {2, 0, 0x4A}, {2, 1, 0x4A}, {2, 2, 0x4A},
  {2, 3, 0x4A}, {2, 4, 0x4A}, {2, 5, 0x4A}
};

void setup() {
  Wire.begin();
  Serial.begin(9600);
  while (!Serial);

  pinMode(S0, OUTPUT);
  pinMode(S1, OUTPUT);
  pinMode(S2, OUTPUT);
  pinMode(S3, OUTPUT);

  Serial.println("Ready for command. Send 'r' to sweep.");
}

void loop() {
  if (Serial.available()) {
    char cmd = Serial.read();
    if (cmd == 'r') {
      runSweep(); // sweep only once
    }
  }
}

void runSweep() {
  Serial.println("Ei,Ej,Vdiff");

  for (int i = 0; i < 22; i++) {
    for (int j = 0; j < 22; j++) {
      if (i == j) continue;

      selectMUX(electrodes[i].channel);
      delayMicroseconds(300);
      int16_t vi = readADS(electrodes[i].adsAddr);

      selectMUX(electrodes[j].channel);
      delayMicroseconds(300);
      int16_t vj = readADS(electrodes[j].adsAddr);

      int16_t vdiff = vi - vj;

      Serial.print(i + 2);
      Serial.print(",");
      Serial.print(j + 2);
      Serial.print(",");
      Serial.println(vdiff);
    }
  }

  Serial.println(" Sweep done. Waiting for next 'r'...");
}

void selectMUX(uint8_t ch) {
  digitalWrite(S0, bitRead(ch, 0));
  digitalWrite(S1, bitRead(ch, 1));
  digitalWrite(S2, bitRead(ch, 2));
  digitalWrite(S3, bitRead(ch, 3));
}

int16_t readADS(uint8_t addr) {
  Wire.beginTransmission(addr);
  Wire.write(0x01);
  Wire.write(0b11000010);
  Wire.write(0b11100011);
  Wire.endTransmission();

  delayMicroseconds(1000);

  Wire.beginTransmission(addr);
  Wire.write(0x00);
  Wire.endTransmission();
  Wire.requestFrom(addr, (uint8_t)2);

  if (Wire.available() < 2) return 0;

  return (Wire.read() << 8) | Wire.read();
}
