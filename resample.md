# Coefficients for oversampling cascades

This is a condensed version of the original "oversampling.xt" included in HIIR.
You only need the "Trans BW" parameters to pass into `coefs_transition`.
See the [iir_4x example](./examples/measure_delay_iir_4x.rs) on how to implement a cascade as specified by the 4x 120 dB below.
Note you always put the largest one on the lowest rate part of the cascade (i.e. on the "outside".)

## 4x

### 4x, 140 dB
```
Total delay : 6.0 spl

Coefficients: 5
Trans BW    : 0.241371

Coefficients: 12
Trans BW    : 0.0330036
```

### 4x, 120 dB
```
Total delay : 5.0 spl

Coefficients: 4
Trans BW    : 0.261666

Coefficients: 10
Trans BW    : 0.0367598
```
### 4x, 99 dB
```
Total delay : 4.0 spl

Coefficients: 3
Trans BW    : 0.283623

Coefficients: 8
Trans BW    : 0.0432267
```
## 8x

### 8x, 133 dB
```
Total delay : 6.0 spl

Coefficients: 3
Trans BW    : 0.373128

Coefficients: 5
Trans BW    : 0.224566

Coefficients: 10
Trans BW    : 0.0507022
```

### 8x, 129 dB
```
Total delay : 6.0 spl

Coefficients: 3
Trans BW    : 0.3901

Coefficients: 4
Trans BW    : 0.284503

Coefficients: 10
Trans BW    : 0.0646551
```

### 8x, 115 dB
```
Total delay : 5.0 spl

Coefficients: 3
Trans BW    : 0.330575

Coefficients: 4
Trans BW    : 0.245813

Coefficients: 8
Trans BW    : 0.0627931
```

### 8x, 96 dB
```
Total delay : 4.0 spl

Coefficients: 2
Trans BW    : 0.384234

Coefficients: 3
Trans BW    : 0.273878

Coefficients: 7
Trans BW    : 0.0544272
```
