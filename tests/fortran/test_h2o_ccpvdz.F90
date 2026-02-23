! Standalone test for libERI: H2O / cc-pVDZ
!
! Exercises s, p, and d kernels (int0000 through int0022 and cross terms).
! Uses the liberi_interface API (no globals, no MDI).
!
! d ordering: x^2, y^2, z^2, xy, xz, yz  (6 Cartesian components)
!
! Usage:
!   ./test_h2o_ccpvdz
!   ./test_h2o_ccpvdz --generate
!
program test_h2o_ccpvdz
  use liberi_types, only: dp, int64
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_destroy, &
                              liberi_setup, liberi_fock_build, liberi_cleanup
  implicit none

  type(liberi_handle_t) :: handle
  integer, parameter :: n_tri = 325
  integer :: i, tri_size
  integer :: nargs, iunit
  character(len=256) :: arg
  logical :: generate
  real(dp) :: max_err, err
  real(dp), parameter :: tol = 1.0d-10
  real(dp), parameter :: pi_val = 3.14159265358979323846264338327950288_dp

  ! Basis set data
  integer :: nsh_val, natoms_val, num_bas_val, mxgtot_val
  integer, allocatable :: ang_mom_arr(:), contr_num_arr(:), sh_loc_arr(:)
  integer, allocatable :: atom_num_arr(:), atom_loc_arr(:)
  integer, allocatable :: start_bas_arr(:), end_bas_arr(:)
  real(dp), allocatable :: exponents_arr(:)
  real(dp), allocatable :: contr_coef_s_arr(:), contr_coef_p_arr(:)
  real(dp), allocatable :: contr_coef_d_arr(:), contr_coef_f_arr(:)
  real(dp), allocatable :: coords_arr(:), schwrz_int_arr(:)
  real(dp), allocatable :: density_arr(:), fock_arr(:)

  ! Reference Fock matrix: H2O/cc-pVDZ, identity density, no screening
  real(dp), parameter :: fa_ref(n_tri) = [ &
                         3.03257788817414493d+01, &  !   1
                         -1.62147486884554626d+00, &  !   2
                         2.09202538203897639d+01, &  !   3
                         1.03668166202817034d+01, &  !   4
                         3.63016223221230447d+01, &  !   5
                         1.86484351856507828d+01, &  !   6
                         0.00000000000000000d+00, &  !   7
                         0.00000000000000000d+00, &  !   8
                         0.00000000000000000d+00, &  !   9
                         2.20128028497883363d+01, &  !  10
                         0.00000000000000000d+00, &  !  11
                         1.38777878078144568d-17, &  !  12
                         0.00000000000000000d+00, &  !  13
                         0.00000000000000000d+00, &  !  14
                         2.21098634796176761d+01, &  !  15
                         7.24023600663820099d-02, &  !  16
                         7.30568658200811516d-01, &  !  17
                         7.51371910495167272d-01, &  !  18
                         0.00000000000000000d+00, &  !  19
                         -1.23599047663347505d-17, &  !  20
                         2.20708723277126460d+01, &  !  21
                         0.00000000000000000d+00, &  !  22
                         0.00000000000000000d+00, &  !  23
                         0.00000000000000000d+00, &  !  24
                         2.95299941637340808d+01, &  !  25
                         0.00000000000000000d+00, &  !  26
                         0.00000000000000000d+00, &  !  27
                         1.62663391345324371d+01, &  !  28
                         0.00000000000000000d+00, &  !  29
                         0.00000000000000000d+00, &  !  30
                         -2.77555756156289135d-17, &  !  31
                         0.00000000000000000d+00, &  !  32
                         2.97713181336327573d+01, &  !  33
                         -3.90312782094781596d-18, &  !  34
                         0.00000000000000000d+00, &  !  35
                         1.64543569439530692d+01, &  !  36
                         -1.91593424735275186d-02, &  !  37
                         5.97601380227355428d-01, &  !  38
                         7.61125967288953209d-01, &  !  39
                         0.00000000000000000d+00, &  !  40
                         -3.46944695195361419d-18, &  !  41
                         2.96743735786140128d+01, &  !  42
                         0.00000000000000000d+00, &  !  43
                         -5.20417042793042128d-18, &  !  44
                         1.63788265153862973d+01, &  !  45
                         5.56318769128527180d+00, &  !  46
                         4.96475386253657689d+01, &  !  47
                         4.53459344149495394d+01, &  !  48
                         0.00000000000000000d+00, &  !  49
                         2.77555756156289135d-17, &  !  50
                         4.23295371715842472d-01, &  !  51
                         0.00000000000000000d+00, &  !  52
                         1.38777878078144568d-17, &  !  53
                         2.91906929249521219d-01, &  !  54
                         6.39852481312498895d+01, &  !  55
                         5.53192106082263280d+00, &  !  56
                         4.99060260801048088d+01, &  !  57
                         4.56106316483041141d+01, &  !  58
                         0.00000000000000000d+00, &  !  59
                         0.00000000000000000d+00, &  !  60
                         5.05060048474446743d-01, &  !  61
                         0.00000000000000000d+00, &  !  62
                         2.22044604925031308d-16, &  !  63
                         4.09110174521912118d-01, &  !  64
                         3.78055527576139596d+01, &  !  65
                         6.45118832262675568d+01, &  !  66
                         5.54448147733827668d+00, &  !  67
                         4.98021866177914490d+01, &  !  68
                         4.55042975965601926d+01, &  !  69
                         0.00000000000000000d+00, &  !  70
                         -5.55111512312578270d-17, &  !  71
                         1.80428776259713208d+00, &  !  72
                         0.00000000000000000d+00, &  !  73
                         0.00000000000000000d+00, &  !  74
                         1.32886562820442933d+00, &  !  75
                         3.77894552314305869d+01, &  !  76
                         3.79362989815979574d+01, &  !  77
                         6.42788773772674631d+01, &  !  78
                         0.00000000000000000d+00, &  !  79
                         0.00000000000000000d+00, &  !  80
                         0.00000000000000000d+00, &  !  81
                         1.11022302462515654d-16, &  !  82
                         0.00000000000000000d+00, &  !  83
                         0.00000000000000000d+00, &  !  84
                         2.77555756156289135d-17, &  !  85
                         0.00000000000000000d+00, &  !  86
                         0.00000000000000000d+00, &  !  87
                         0.00000000000000000d+00, &  !  88
                         0.00000000000000000d+00, &  !  89
                         0.00000000000000000d+00, &  !  90
                         6.68879400295607809d+01, &  !  91
                         0.00000000000000000d+00, &  !  92
                         0.00000000000000000d+00, &  !  93
                         0.00000000000000000d+00, &  !  94
                         1.15361008198086123d+00, &  !  95
                         0.00000000000000000d+00, &  !  96
                         0.00000000000000000d+00, &  !  97
                         8.37306518645006714d-01, &  !  98
                         0.00000000000000000d+00, &  !  99
                         0.00000000000000000d+00, &  ! 100
                         0.00000000000000000d+00, &  ! 101
                         0.00000000000000000d+00, &  ! 102
                         0.00000000000000000d+00, &  ! 103
                         0.00000000000000000d+00, &  ! 104
                         6.67682283149877520d+01, &  ! 105
                         0.00000000000000000d+00, &  ! 106
                         4.16333634234433703d-17, &  ! 107
                         -1.38777878078144568d-17, &  ! 108
                         0.00000000000000000d+00, &  ! 109
                         1.29523065639120971d+00, &  ! 110
                         5.55111512312578270d-17, &  ! 111
                         0.00000000000000000d+00, &  ! 112
                         1.04030849426874483d+00, &  ! 113
                         -2.77555756156289135d-17, &  ! 114
                         0.00000000000000000d+00, &  ! 115
                         0.00000000000000000d+00, &  ! 116
                         0.00000000000000000d+00, &  ! 117
                         0.00000000000000000d+00, &  ! 118
                         0.00000000000000000d+00, &  ! 119
                         6.72263848441013607d+01, &  ! 120
                         1.59739777664195870d+00, &  ! 121
                         9.22765098009861262d+00, &  ! 122
                         1.03791766666979459d+01, &  ! 123
                         0.00000000000000000d+00, &  ! 124
                         8.25533100097732131d+00, &  ! 125
                         6.49814794642948090d+00, &  ! 126
                         0.00000000000000000d+00, &  ! 127
                         9.71140033564597438d+00, &  ! 128
                         7.62991534722465392d+00, &  ! 129
                         7.87828590066039247d+00, &  ! 130
                         1.39123181440287400d+01, &  ! 131
                         1.16419009904662101d+01, &  ! 132
                         0.00000000000000000d+00, &  ! 133
                         0.00000000000000000d+00, &  ! 134
                         8.10469220196533691d+00, &  ! 135
                         5.81249650039644727d+00, &  ! 136
                         3.55744472490763197d+00, &  ! 137
                         1.90017552730139734d+01, &  ! 138
                         2.18099418624883228d+01, &  ! 139
                         0.00000000000000000d+00, &  ! 140
                         6.46437925252116319d+00, &  ! 141
                         5.42182629075120559d+00, &  ! 142
                         0.00000000000000000d+00, &  ! 143
                         9.00013035679119255d+00, &  ! 144
                         7.51288791564411174d+00, &  ! 145
                         2.09992055257035766d+01, &  ! 146
                         2.21581613923269316d+01, &  ! 147
                         2.18121148243743725d+01, &  ! 148
                         0.00000000000000000d+00, &  ! 149
                         0.00000000000000000d+00, &  ! 150
                         1.51316421968093939d+00, &  ! 151
                         1.23759750747926383d+01, &  ! 152
                         1.30362925903957425d+01, &  ! 153
                         0.00000000000000000d+00, &  ! 154
                         0.00000000000000000d+00, &  ! 155
                         0.00000000000000000d+00, &  ! 156
                         1.17157710703045268d+01, &  ! 157
                         0.00000000000000000d+00, &  ! 158
                         0.00000000000000000d+00, &  ! 159
                         1.32107533889090742d+01, &  ! 160
                         0.00000000000000000d+00, &  ! 161
                         0.00000000000000000d+00, &  ! 162
                         0.00000000000000000d+00, &  ! 163
                         0.00000000000000000d+00, &  ! 164
                         0.00000000000000000d+00, &  ! 165
                         1.73551318114398008d+01, &  ! 166
                         1.36639294407713603d+01, &  ! 167
                         0.00000000000000000d+00, &  ! 168
                         0.00000000000000000d+00, &  ! 169
                         0.00000000000000000d+00, &  ! 170
                         1.63368062182103628d+01, &  ! 171
                         -4.69306421067425283d+00, &  ! 172
                         -1.47088327360470092d+01, &  ! 173
                         -1.33721451122969484d+01, &  ! 174
                         0.00000000000000000d+00, &  ! 175
                         -5.50923395265970317d+00, &  ! 176
                         -1.35566001951841848d+01, &  ! 177
                         0.00000000000000000d+00, &  ! 178
                         6.25805693030418719d-02, &  ! 179
                         -1.03445229844505970d+01, &  ! 180
                         -1.65805707156567017d+01, &  ! 181
                         -1.59728001288368020d+01, &  ! 182
                         -2.88727977121916624d+01, &  ! 183
                         0.00000000000000000d+00, &  ! 184
                         0.00000000000000000d+00, &  ! 185
                         -1.28558595255611099d+01, &  ! 186
                         -2.53941953582625679d+00, &  ! 187
                         -3.43775930401710150d+00, &  ! 188
                         0.00000000000000000d+00, &  ! 189
                         1.70394817950801247d+01, &  ! 190
                         -3.63760316848446275d+00, &  ! 191
                         -1.11839123435888546d+01, &  ! 192
                         -1.00994513345176280d+01, &  ! 193
                         0.00000000000000000d+00, &  ! 194
                         -1.31522192053660394d+01, &  ! 195
                         1.42248048099394131d+00, &  ! 196
                         0.00000000000000000d+00, &  ! 197
                         -9.93087388377740865d+00, &  ! 198
                         5.46511269699949231d+00, &  ! 199
                         -1.27632905381695210d+01, &  ! 200
                         -2.79388795832668535d+01, &  ! 201
                         -6.12861635712836961d+00, &  ! 202
                         0.00000000000000000d+00, &  ! 203
                         0.00000000000000000d+00, &  ! 204
                         -2.93853011868754344d+00, &  ! 205
                         -1.73785614100221841d+00, &  ! 206
                         -2.33671161772376124d+00, &  ! 207
                         0.00000000000000000d+00, &  ! 208
                         8.97646127469533428d-01, &  ! 209
                         1.67170318126069226d+01, &  ! 210
                         1.59739777664195870d+00, &  ! 211
                         9.22765098009861262d+00, &  ! 212
                         1.03791766666979459d+01, &  ! 213
                         0.00000000000000000d+00, &  ! 214
                         -8.25533100097731953d+00, &  ! 215
                         6.49814794642948090d+00, &  ! 216
                         0.00000000000000000d+00, &  ! 217
                         -9.71140033564597438d+00, &  ! 218
                         7.62991534722465481d+00, &  ! 219
                         7.87828590066039247d+00, &  ! 220
                         1.39123181440287400d+01, &  ! 221
                         1.16419009904662119d+01, &  ! 222
                         0.00000000000000000d+00, &  ! 223
                         0.00000000000000000d+00, &  ! 224
                         -8.10469220196533691d+00, &  ! 225
                         1.42833685344968142d+00, &  ! 226
                         5.84532353761637147d+00, &  ! 227
                         0.00000000000000000d+00, &  ! 228
                         -3.50274435662834627d+00, &  ! 229
                         -1.42141328672058059d-01, &  ! 230
                         5.81249650039644816d+00, &  ! 231
                         3.55744472490763197d+00, &  ! 232
                         1.90017552730139734d+01, &  ! 233
                         2.18099418624883228d+01, &  ! 234
                         0.00000000000000000d+00, &  ! 235
                         -6.46437925252116319d+00, &  ! 236
                         5.42182629075120559d+00, &  ! 237
                         0.00000000000000000d+00, &  ! 238
                         -9.00013035679119078d+00, &  ! 239
                         7.51288791564411262d+00, &  ! 240
                         2.09992055257035766d+01, &  ! 241
                         2.21581613923269316d+01, &  ! 242
                         2.18121148243743725d+01, &  ! 243
                         0.00000000000000000d+00, &  ! 244
                         0.00000000000000000d+00, &  ! 245
                         -1.51316421968093939d+00, &  ! 246
                         5.84532353761637236d+00, &  ! 247
                         1.70917531953688595d+01, &  ! 248
                         0.00000000000000000d+00, &  ! 249
                         -7.05724770764219933d+00, &  ! 250
                         -1.07455201460047145d+00, &  ! 251
                         1.23759750747926383d+01, &  ! 252
                         1.30362925903957425d+01, &  ! 253
                         0.00000000000000000d+00, &  ! 254
                         0.00000000000000000d+00, &  ! 255
                         0.00000000000000000d+00, &  ! 256
                         1.17157710703045268d+01, &  ! 257
                         0.00000000000000000d+00, &  ! 258
                         0.00000000000000000d+00, &  ! 259
                         1.32107533889090742d+01, &  ! 260
                         0.00000000000000000d+00, &  ! 261
                         0.00000000000000000d+00, &  ! 262
                         0.00000000000000000d+00, &  ! 263
                         0.00000000000000000d+00, &  ! 264
                         0.00000000000000000d+00, &  ! 265
                         -1.73551318114398008d+01, &  ! 266
                         1.36639294407713603d+01, &  ! 267
                         0.00000000000000000d+00, &  ! 268
                         0.00000000000000000d+00, &  ! 269
                         0.00000000000000000d+00, &  ! 270
                         1.82826880160592231d+00, &  ! 271
                         0.00000000000000000d+00, &  ! 272
                         0.00000000000000000d+00, &  ! 273
                         0.00000000000000000d+00, &  ! 274
                         0.00000000000000000d+00, &  ! 275
                         1.63368062182103628d+01, &  ! 276
                         4.69306421067425372d+00, &  ! 277
                         1.47088327360470092d+01, &  ! 278
                         1.33721451122969484d+01, &  ! 279
                         0.00000000000000000d+00, &  ! 280
                         -5.50923395265970228d+00, &  ! 281
                         1.35566001951841866d+01, &  ! 282
                         0.00000000000000000d+00, &  ! 283
                         6.25805693030419830d-02, &  ! 284
                         1.03445229844505988d+01, &  ! 285
                         1.65805707156567017d+01, &  ! 286
                         1.59728001288368020d+01, &  ! 287
                         2.88727977121916624d+01, &  ! 288
                         0.00000000000000000d+00, &  ! 289
                         0.00000000000000000d+00, &  ! 290
                         -1.28558595255611099d+01, &  ! 291
                         3.50274435662834627d+00, &  ! 292
                         7.05724770764219933d+00, &  ! 293
                         0.00000000000000000d+00, &  ! 294
                         -1.01118436154760367d+01, &  ! 295
                         -3.99328952664650794d-01, &  ! 296
                         2.53941953582625679d+00, &  ! 297
                         3.43775930401710239d+00, &  ! 298
                         0.00000000000000000d+00, &  ! 299
                         1.70394817950801247d+01, &  ! 300
                         -3.63760316848446275d+00, &  ! 301
                         -1.11839123435888563d+01, &  ! 302
                         -1.00994513345176298d+01, &  ! 303
                         0.00000000000000000d+00, &  ! 304
                         1.31522192053660429d+01, &  ! 305
                         1.42248048099394087d+00, &  ! 306
                         0.00000000000000000d+00, &  ! 307
                         9.93087388377741043d+00, &  ! 308
                         5.46511269699949320d+00, &  ! 309
                         -1.27632905381695210d+01, &  ! 310
                         -2.79388795832668499d+01, &  ! 311
                         -6.12861635712836872d+00, &  ! 312
                         0.00000000000000000d+00, &  ! 313
                         0.00000000000000000d+00, &  ! 314
                         2.93853011868754344d+00, &  ! 315
                         -1.42141328672058059d-01, &  ! 316
                         -1.07455201460047145d+00, &  ! 317
                         0.00000000000000000d+00, &  ! 318
                         3.99328952664650627d-01, &  ! 319
                         1.76256767846435514d+00, &  ! 320
                         -1.73785614100221841d+00, &  ! 321
                         -2.33671161772376212d+00, &  ! 322
                         0.00000000000000000d+00, &  ! 323
                         -8.97646127469533428d-01, &  ! 324
                         1.67170318126069226d+01 &  ! 325
                         ]

  ! cc-pVDZ basis for oxygen (BSE, general contraction):
  !   S  9  2.00   (general contraction → 2 contracted s functions)
  !       11720.0    0.00071   -0.00016
  !        1759.0    0.00547   -0.001263
  !         400.8    0.027837  -0.006267
  !         113.7    0.1048    -0.025716
  !          37.03   0.283062  -0.070924
  !          13.27   0.448719  -0.165411
  !           5.025  0.270952  -0.116955
  !           1.013  0.015458   0.557368
  !           0.3023 -0.002585  0.572759
  !   S  1  1.00
  !           0.3023  1.0
  !   P  4  1.00   (general contraction)
  !          17.70   0.043018
  !           3.854  0.228913
  !           1.046  0.508728
  !           0.2753 0.460531
  !   P  1  1.00
  !           0.2753  1.0
  !   D  1  1.00
  !           1.185   1.0
  !
  ! cc-pVDZ basis for hydrogen:
  !   S  3  1.00
  !          13.01    0.019685
  !           1.962   0.137977
  !           0.4446  0.478148
  !   S  1  1.00
  !           0.122   1.0
  !   P  1  1.00
  !           0.727   1.0

  ! Parse --generate flag
  generate = .false.
  nargs = command_argument_count()
  do i = 1, nargs
    call get_command_argument(i, arg)
    if (trim(arg) == '--generate') generate = .true.
  end do

  ! ============================================
  ! Molecule: H2O
  !   O at origin, H atoms symmetric in yz-plane
  !   r(OH) = 1.8096 bohr, angle(HOH) = 104.51 deg
  ! ============================================
  natoms_val = 3
  nsh_val = 12         ! 6 shells on O, 3 on each H
  num_bas_val = 25     ! O: 1+1+1+3+3+6=15, H: 1+1+3=5, total 15+5+5=25
  mxgtot_val = 9       ! max primitives in any shell (O general-contracted s)
  tri_size = num_bas_val * (num_bas_val + 1) / 2   ! 325

  ! ============================================
  ! Allocate arrays
  ! ============================================
  allocate(ang_mom_arr(nsh_val))
  allocate(contr_num_arr(nsh_val))
  allocate(sh_loc_arr(nsh_val))
  allocate(atom_num_arr(nsh_val))
  allocate(atom_loc_arr(nsh_val))
  allocate(start_bas_arr(nsh_val))
  allocate(end_bas_arr(nsh_val))
  allocate(coords_arr(3 * natoms_val))
  allocate(exponents_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_s_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_p_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_d_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_f_arr(mxgtot_val * nsh_val))
  allocate(schwrz_int_arr(nsh_val * (nsh_val + 1) / 2))
  allocate(density_arr(tri_size))
  allocate(fock_arr(tri_size))

  ! Initialize coefficient arrays
  exponents_arr = 0.0_dp
  contr_coef_s_arr = 0.0_dp
  contr_coef_p_arr = 0.0_dp
  contr_coef_d_arr = 0.0_dp
  contr_coef_f_arr = 0.0_dp

  ! ============================================
  ! Atomic coordinates (bohr)
  !   O at origin, H1 and H2 in yz-plane
  ! ============================================
  coords_arr(1) = 0.0_dp; coords_arr(2) = 0.0_dp; coords_arr(3) = 0.0_dp
  coords_arr(4) = 0.0_dp; coords_arr(5) = 1.431390_dp; coords_arr(6) = 1.107160_dp
  coords_arr(7) = 0.0_dp; coords_arr(8) = -1.431390_dp; coords_arr(9) = 1.107160_dp

  ! ============================================
  ! Basis set: cc-pVDZ for O (atom 1)
  ! ============================================

  ! Shell 1: s, 9 primitives (general contraction, 1st function)
  ang_mom_arr(1) = 1
  contr_num_arr(1) = 9
  sh_loc_arr(1) = 1
  atom_num_arr(1) = 1
  atom_loc_arr(1) = 1
  start_bas_arr(1) = 1; end_bas_arr(1) = 1
  exponents_arr(1) = 11720.0000000_dp
  exponents_arr(2) = 1759.0000000_dp
  exponents_arr(3) = 400.8000000_dp
  exponents_arr(4) = 113.7000000_dp
  exponents_arr(5) = 37.0300000_dp
  exponents_arr(6) = 13.2700000_dp
  exponents_arr(7) = 5.0250000_dp
  exponents_arr(8) = 1.0130000_dp
  exponents_arr(9) = 0.3023000_dp
  contr_coef_s_arr(1) = 0.0007100_dp*gto_norm_s(11720.0000000_dp)
  contr_coef_s_arr(2) = 0.0054700_dp*gto_norm_s(1759.0000000_dp)
  contr_coef_s_arr(3) = 0.0278370_dp*gto_norm_s(400.8000000_dp)
  contr_coef_s_arr(4) = 0.1048000_dp*gto_norm_s(113.7000000_dp)
  contr_coef_s_arr(5) = 0.2830620_dp*gto_norm_s(37.0300000_dp)
  contr_coef_s_arr(6) = 0.4487190_dp*gto_norm_s(13.2700000_dp)
  contr_coef_s_arr(7) = 0.2709520_dp*gto_norm_s(5.0250000_dp)
  contr_coef_s_arr(8) = 0.0154580_dp*gto_norm_s(1.0130000_dp)
  contr_coef_s_arr(9) = -0.0025850_dp*gto_norm_s(0.3023000_dp)

  ! Shell 2: s, 9 primitives (general contraction, 2nd function)
  ang_mom_arr(2) = 1
  contr_num_arr(2) = 9
  sh_loc_arr(2) = 10
  atom_num_arr(2) = 1
  atom_loc_arr(2) = 2
  start_bas_arr(2) = 2; end_bas_arr(2) = 2
  exponents_arr(10) = 11720.0000000_dp
  exponents_arr(11) = 1759.0000000_dp
  exponents_arr(12) = 400.8000000_dp
  exponents_arr(13) = 113.7000000_dp
  exponents_arr(14) = 37.0300000_dp
  exponents_arr(15) = 13.2700000_dp
  exponents_arr(16) = 5.0250000_dp
  exponents_arr(17) = 1.0130000_dp
  exponents_arr(18) = 0.3023000_dp
  contr_coef_s_arr(10) = -0.0001600_dp*gto_norm_s(11720.0000000_dp)
  contr_coef_s_arr(11) = -0.0012630_dp*gto_norm_s(1759.0000000_dp)
  contr_coef_s_arr(12) = -0.0062670_dp*gto_norm_s(400.8000000_dp)
  contr_coef_s_arr(13) = -0.0257160_dp*gto_norm_s(113.7000000_dp)
  contr_coef_s_arr(14) = -0.0709240_dp*gto_norm_s(37.0300000_dp)
  contr_coef_s_arr(15) = -0.1654110_dp*gto_norm_s(13.2700000_dp)
  contr_coef_s_arr(16) = -0.1169550_dp*gto_norm_s(5.0250000_dp)
  contr_coef_s_arr(17) = 0.5573680_dp*gto_norm_s(1.0130000_dp)
  contr_coef_s_arr(18) = 0.5727590_dp*gto_norm_s(0.3023000_dp)

  ! Shell 3: s, 1 primitive (uncontracted)
  ang_mom_arr(3) = 1
  contr_num_arr(3) = 1
  sh_loc_arr(3) = 19
  atom_num_arr(3) = 1
  atom_loc_arr(3) = 3
  start_bas_arr(3) = 3; end_bas_arr(3) = 3
  exponents_arr(19) = 0.3023000_dp
  contr_coef_s_arr(19) = 1.0000000_dp*gto_norm_s(0.3023000_dp)

  ! Shell 4: p, 4 primitives (general contraction)
  ang_mom_arr(4) = 2
  contr_num_arr(4) = 4
  sh_loc_arr(4) = 20
  atom_num_arr(4) = 1
  atom_loc_arr(4) = 4
  start_bas_arr(4) = 4; end_bas_arr(4) = 6
  exponents_arr(20) = 17.7000000_dp
  exponents_arr(21) = 3.8540000_dp
  exponents_arr(22) = 1.0460000_dp
  exponents_arr(23) = 0.2753000_dp
  contr_coef_p_arr(20) = 0.0430180_dp*gto_norm_p(17.7000000_dp)
  contr_coef_p_arr(21) = 0.2289130_dp*gto_norm_p(3.8540000_dp)
  contr_coef_p_arr(22) = 0.5087280_dp*gto_norm_p(1.0460000_dp)
  contr_coef_p_arr(23) = 0.4605310_dp*gto_norm_p(0.2753000_dp)

  ! Shell 5: p, 1 primitive (uncontracted)
  ang_mom_arr(5) = 2
  contr_num_arr(5) = 1
  sh_loc_arr(5) = 24
  atom_num_arr(5) = 1
  atom_loc_arr(5) = 7
  start_bas_arr(5) = 7; end_bas_arr(5) = 9
  exponents_arr(24) = 0.2753000_dp
  contr_coef_p_arr(24) = 1.0000000_dp*gto_norm_p(0.2753000_dp)

  ! Shell 6: d, 1 primitive (uncontracted, 6 Cartesian components)
  ang_mom_arr(6) = 3
  contr_num_arr(6) = 1
  sh_loc_arr(6) = 25
  atom_num_arr(6) = 1
  atom_loc_arr(6) = 10
  start_bas_arr(6) = 10; end_bas_arr(6) = 15
  exponents_arr(25) = 1.1850000_dp
  contr_coef_d_arr(25) = 1.0000000_dp*gto_norm_d(1.1850000_dp)

  ! ============================================
  ! Basis set: cc-pVDZ for H (atom 2)
  ! ============================================

  ! Shell 7: s, 3 primitives
  ang_mom_arr(7) = 1
  contr_num_arr(7) = 3
  sh_loc_arr(7) = 26
  atom_num_arr(7) = 2
  atom_loc_arr(7) = 16
  start_bas_arr(7) = 16; end_bas_arr(7) = 16
  exponents_arr(26) = 13.0100000_dp
  exponents_arr(27) = 1.9620000_dp
  exponents_arr(28) = 0.4446000_dp
  contr_coef_s_arr(26) = 0.0196850_dp*gto_norm_s(13.0100000_dp)
  contr_coef_s_arr(27) = 0.1379770_dp*gto_norm_s(1.9620000_dp)
  contr_coef_s_arr(28) = 0.4781480_dp*gto_norm_s(0.4446000_dp)

  ! Shell 8: s, 1 primitive
  ang_mom_arr(8) = 1
  contr_num_arr(8) = 1
  sh_loc_arr(8) = 29
  atom_num_arr(8) = 2
  atom_loc_arr(8) = 17
  start_bas_arr(8) = 17; end_bas_arr(8) = 17
  exponents_arr(29) = 0.1220000_dp
  contr_coef_s_arr(29) = 1.0000000_dp*gto_norm_s(0.1220000_dp)

  ! Shell 9: p, 1 primitive
  ang_mom_arr(9) = 2
  contr_num_arr(9) = 1
  sh_loc_arr(9) = 30
  atom_num_arr(9) = 2
  atom_loc_arr(9) = 18
  start_bas_arr(9) = 18; end_bas_arr(9) = 20
  exponents_arr(30) = 0.7270000_dp
  contr_coef_p_arr(30) = 1.0000000_dp*gto_norm_p(0.7270000_dp)

  ! ============================================
  ! Basis set: cc-pVDZ for H (atom 3)
  ! ============================================

  ! Shell 10: s, 3 primitives
  ang_mom_arr(10) = 1
  contr_num_arr(10) = 3
  sh_loc_arr(10) = 31
  atom_num_arr(10) = 3
  atom_loc_arr(10) = 21
  start_bas_arr(10) = 21; end_bas_arr(10) = 21
  exponents_arr(31) = 13.0100000_dp
  exponents_arr(32) = 1.9620000_dp
  exponents_arr(33) = 0.4446000_dp
  contr_coef_s_arr(31) = 0.0196850_dp*gto_norm_s(13.0100000_dp)
  contr_coef_s_arr(32) = 0.1379770_dp*gto_norm_s(1.9620000_dp)
  contr_coef_s_arr(33) = 0.4781480_dp*gto_norm_s(0.4446000_dp)

  ! Shell 11: s, 1 primitive
  ang_mom_arr(11) = 1
  contr_num_arr(11) = 1
  sh_loc_arr(11) = 34
  atom_num_arr(11) = 3
  atom_loc_arr(11) = 22
  start_bas_arr(11) = 22; end_bas_arr(11) = 22
  exponents_arr(34) = 0.1220000_dp
  contr_coef_s_arr(34) = 1.0000000_dp*gto_norm_s(0.1220000_dp)

  ! Shell 12: p, 1 primitive
  ang_mom_arr(12) = 2
  contr_num_arr(12) = 1
  sh_loc_arr(12) = 35
  atom_num_arr(12) = 3
  atom_loc_arr(12) = 23
  start_bas_arr(12) = 23; end_bas_arr(12) = 25
  exponents_arr(35) = 0.7270000_dp
  contr_coef_p_arr(35) = 1.0000000_dp*gto_norm_p(0.7270000_dp)

  ! ============================================
  ! Disable Schwarz screening (set large values)
  ! ============================================
  schwrz_int_arr = 1.0d10

  ! ============================================
  ! Density matrix — identity (diagonal = 1)
  ! Triangular packed index: diag(k) = k*(k+1)/2
  ! ============================================
  density_arr = 0.0_dp
  do i = 1, num_bas_val
    density_arr(i * (i + 1) / 2) = 1.0_dp
  end do
  fock_arr = 0.0_dp

  ! ============================================
  ! Use liberi_interface API
  ! ============================================
  print *, "Creating handle..."
  call liberi_create(handle)

  print *, "Setting up basis and computing shell pairs..."
  call liberi_setup(handle, nsh_val, natoms_val, num_bas_val, mxgtot_val, &
                    ang_mom_arr, contr_num_arr, sh_loc_arr,               &
                    atom_num_arr, atom_loc_arr,                           &
                    start_bas_arr, end_bas_arr,                           &
                    exponents_arr,                                        &
                    contr_coef_s_arr, contr_coef_p_arr,                   &
                    contr_coef_d_arr, contr_coef_f_arr,                   &
                    coords_arr, schwrz_int_arr,                           &
                    0, 1)  ! my_rank=0, num_procs=1

  print *, "Computing Fock matrix..."
  call liberi_fock_build(handle, density_arr, fock_arr, tri_size)

  ! ============================================
  ! Generate or validate Fock matrix
  ! ============================================
  if (generate) then
    iunit = 20
    open(unit=iunit, file='fa_h2o_reference.dat', form='formatted', status='replace')
    write(iunit, '(I6)') tri_size
    do i = 1, tri_size
      write(iunit, '(ES26.17)') fock_arr(i)
    end do
    close(iunit)
    print *, ""
    print *, "Wrote reference to: fa_h2o_reference.dat"
    do i = 1, tri_size
      print '(A,I4,A,ES24.15)', "  fa(", i, ") = ", fock_arr(i)
    end do
  else
    max_err = 0.0_dp
    do i = 1, tri_size
      err = abs(fock_arr(i) - fa_ref(i))
      if (err > max_err) max_err = err
      if (err > tol) then
        print '(A,I4,A,ES24.15,A,ES24.15,A,ES10.2)', &
          "  FAIL fa(", i, ") = ", fock_arr(i), "  ref = ", fa_ref(i), "  err = ", err
      end if
    end do

    print *, ""
    if (max_err <= tol) then
      print '(A,ES10.2)', " PASS - max error: ", max_err
    else
      print '(A,ES10.2)', " FAIL - max error: ", max_err
      call liberi_destroy(handle)
      error stop 1
    end if
  end if

  ! ============================================
  ! Cleanup
  ! ============================================
  print *, "Cleaning up..."
  call liberi_cleanup(handle)
  call liberi_destroy(handle)

  deallocate(ang_mom_arr, contr_num_arr, sh_loc_arr)
  deallocate(atom_num_arr, atom_loc_arr, start_bas_arr, end_bas_arr)
  deallocate(coords_arr, exponents_arr)
  deallocate(contr_coef_s_arr, contr_coef_p_arr, contr_coef_d_arr, contr_coef_f_arr)
  deallocate(schwrz_int_arr, density_arr, fock_arr)

  print *, "Done!"

contains

  ! GAMESS normalization for s-type primitive: N = (2α/π)^(3/4)
  pure function gto_norm_s(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = (2.0_dp*alpha/pi_val)**0.75_dp
  end function

  ! GAMESS normalization for p-type primitive: N = 2√α × (2α/π)^(3/4)
  pure function gto_norm_p(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 2.0_dp*sqrt(alpha)*(2.0_dp*alpha/pi_val)**0.75_dp
  end function

  ! GAMESS normalization for d-type primitive: N = 4α × (2α/π)^(3/4)
  ! (normalized for off-diagonal components xy, xz, yz;
  !  diagonal components xx, yy, zz have extra 1/√3 handled in kernel)
  pure function gto_norm_d(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 4.0_dp*alpha*(2.0_dp*alpha/pi_val)**0.75_dp
  end function

end program test_h2o_ccpvdz
