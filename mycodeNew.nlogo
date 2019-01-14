;Robert Klock and Carlo Maley
;Autumn 2016 to Winter 2017
;Last updated Spring 2017
;A model of double-bind therapy

globals
[chemo_strength
antibody_strength
die?
both_inactive
max_pop
anti births
population
max_blue
max_green
cc_patches
;antibodypop
immune_predation_deactivated
immune_predation_activated
green-count
blue-count
pd1pop
ctla4pop
pd1-pop
ctla4-pop
justpd1
justctla4
mutationrate
mutationsize
bothpop
averagesurvivalrate
averageproliferationrate
]

turtles-own[
   rate_proliferation
   rate_survival
   CTLA4
   survival_probability
   PD1]

to setup
  clear-all
  reset-ticks
  set cc_patches 0 ;cc patches keeps track of the total number of patches in the simulation at carrying capacitymn,
  set mutationrate .05
  set mutationsize 10
  set max_pop 0
  set both_inactive 0
  set max_blue 0
  set max_green 0
  set immune_predation_deactivated 55
  set immune_predation_activated 10
  set population count patches * 0.90
  ask n-of population patches [sprout 1 [set shape "circle 2"
  set size 0.8
  setxy pxcor pycor
  let var random-float 1
  ifelse( var < 0.5)[
      set rate_proliferation random 100
      set rate_survival 110 - rate_proliferation
    ]
    [
      set rate_survival random 100
      set rate_proliferation 110 - rate_survival
    ]
  ifelse rate_proliferation > rate_survival
          [
            set color scale-color green rate_proliferation 170 0
          ]
          [
            set color scale-color blue rate_survival 170 0
         ]
      ]
    ]

  ask patches [
    if turtles-here = cc [set cc_patches cc_patches + 1]
  ]
ask turtles [
  ifelse
  (random-float 1 < mutationrate) [set pd1 1] [set pd1 0]
  ]


ask turtles [
  ifelse
  (random-float 1 < mutationrate) [set ctla4 1] [set ctla4 0]
  ]

ask turtles [set survival_probability (rate_survival * 0.01) ] ;if rate survival is 77, survival prob is 0.77

;ask (n-of (%t * population) turtles)[set antibodies_ppresent 1]
end

to update
  set green-count count turtles with [rate_proliferation > 50]
  set blue-count count turtles with [rate_survival > 50]
end

to setupPop
  ask turtles
   [
    proliferate
    reach_capacity
    ]
end

to checkDeath
if Ctla4 = 0 and pd1 = 0
  [
   if random-float 1 < immune_susceptibility [die]
  ]
let a random-float 1
let apoptosis random 100
ifelse apoptosis < apoptosis_rate  = true [die]
[
;A CTLA4 mutant (resistant to PD1 antibody)
;cell dies due to:
;Background apoptosis rate/ Chemotherapy - it must be super-sensitive to chemo in order to set up the double-bind

;if pd1 is activated
; Immune_susceptibility -> background chance that a cancer cell is killed by the immune system

if immuno_switch
[ ;originally had the ctla4/pd1 code above
  if pd1 = 1[
   if random-float 1 < immune_susceptibility + ((1 - immune_susceptibility) * survival_probability * doublebind_immune)[die]
  ]

]
;CHEMO PREDATION
;resistance to one = sensitivity to other
; The bug (2/3/17) has to do with the amount of immune susceptibility
if chemo_switch [
  let b random-float 1
  if  CTLA4 = 1 [
    if b < (1 - survival_probability + (doublebind_chemo * survival_probability)) * (chemo_dose) [die]
    ;if b  > (1 - chemo_dose) + (survival_probability * (1 - doublebind_chemo)) [die]
  ]
  if CTLA4 = 0 [
    if b < (1 - survival_probability) * (chemo_dose) [die]
  ]
  ]
]
end
;selection check experiments
;might have to change chemotherapy for bigger population selectivity
;average rate of proliferation (t-test plot) when there is no chemo therapy the average rate of proliferation is x and when chemo is on the average proliferation is hopefully lower
;when there is chemo, rate should be lower (one-sided t test)
;immmune predation test : when immune susceptibility is low/high, the amount of cells that are Pd1 or ctla4 activated should be neutral (half) /high respectively
;
;If CTLA4 [
;    If random number > sensitivity_factor * survival_probability, then die
;  ]
;  Else [
;    If random number > survival_probability, then die

to go
;  ask turtles [proliferate]
;  if population = 800
;  [
  set population count turtles
  if population > max_pop [set max_pop population]
  ask turtles
  [;if random 100 < apoptosis_rate [die]
    ifelse count turtles-here < cc  ; basically coin flip to tell if cell is going to proliferate or not
    [ if rate_proliferation + (max_pro * (100 - rate_proliferation)) > random 100 [
        ;originally "if rate_proliferation * max_pro > random 100...
      proliferate ]]
    [reach_capacity
     ]

  ]
  ask turtles
  [
    checkdeath
    update-plot
  ]
  if timed_therapy[
  if ticks >= chemo_start [set chemo_switch true]
  if ticks >= immuno_start [set immuno_switch true]
  ]
  if population_therapy[
    if population > population_start
    [if chemo_start > immuno_start[
      if ticks <= chemo_start[
        if ticks >= immuno_start [
          set immuno_switch true]]
      if ticks > chemo_start[
        set chemo_switch true
        set immuno_switch false]
      ]

      if immuno_start > chemo_start[
        if ticks <= immuno_start[
          if ticks >= chemo_start[
           set  chemo_switch true]]
        if ticks > immuno_start[
          set chemo_switch false
          set immuno_switch true]
  ]
  ]
  ]


  if ticks = 2400 [stop]
  if count turtles > 0 [
    set averagesurvivalrate mean ([survival_probability] of turtles)
    set averageproliferationrate mean ([rate_proliferation] of turtles)
  ]
  tick

end

to update-plot

  ;set-current-plot _clarify-duplicate-plot-name "population"
  set green-count count turtles with [rate_proliferation > rate_survival]
  set blue-count count turtles with [rate_survival > rate_proliferation]
  set bothpop count turtles with [pd1 = 1 and ctla4 = 1]
  set pd1pop count turtles with [PD1 = 1]
  set ctla4pop count turtles with [CTLA4 = 1]
  set justpd1 count turtles with [ctla4 = 0 and pd1 = 1]
  set justctla4 count turtles with [ctla4 = 1 and pd1 = 0]
  set pd1-pop count turtles with [pd1 = 0]
  set ctla4-pop count turtles with [ctla4 = 0]
  set both_inactive count turtles with [CTLA4 = 0 and PD1 = 0]


end

to reach_capacity ; if cells have a neighbor patch that has not reached cc, cells can move to empty space
      let num_turtles count turtles-here with [self != myself] ; count number of turtles on patch
      let potential_spots neighbors
      ;let potential_spots neighbors4 ;moves turtles in just 4 directions - N, E, S, W
      let chance count potential_spots with [count turtles-here < cc] ;see if neighbors have empty spot
      if chance > 0 [
        let target_spot one-of potential_spots with [count turtles-here < cc] ; if so, move to empty spot
        face target_spot
        fd 1
      ]


end


to proliferate
  hatch 1
    [
      set shape "circle 2"
      set size 0.8

      ;set a_rate apoptosis_rate
      ;set ticks-since-here 0 ; set down counter for time to proliferate again. if proliferated - must wait 24 ticks before proliferating again.

      ;set rate_proliferation random 100
      if random-float 1 < mutationrate [
        let change random mutationsize
        set change change - ( .5 * mutationsize)
        ;set rate_proliferation rate_proliferation + change
        set rate_survival rate_survival + change
        let var random-float 1
  ifelse( var < 0.5)[
      set rate_proliferation random 100
      set rate_survival 110 - rate_proliferation
    ]
    [
      set rate_survival random 100
      set rate_proliferation 110 - rate_survival
    ]
        if (rate_proliferation > 100)[set rate_proliferation 99]
        if (rate_survival > 100) [set rate_survival 99]
  ifelse rate_proliferation > rate_survival
          [
            set color scale-color green rate_proliferation 170 0
          ]
          [
            set color scale-color blue rate_survival 170 0
         ]
      ]
      ;set rate_proliferation 100 - rate_survival
     ; set rate_survival 100 - rate_proliferation
      set survival_probability (rate_survival * 0.01)

      ifelse rate_proliferation > rate_survival
      [
        set color scale-color green rate_proliferation 170 0

      ]
      [
        set color scale-color blue rate_survival 170 0

      ]
      ;The following code didn't do anything, we set ctla and pd1 randomly afterwards anyays
      if random-float 1 < mutationrate [
        ifelse CTLA4 = 1
        [set CTLA4 0]
        [set CTLA4 1]
      ]
      if random-float 1 < mutationrate [
        ifelse PD1 = 1
        [set PD1 0]
        [set PD1 1]
      ]

     ; if pd1 = 0 [
     ; set chemo_sensitivity abs (( chemo_sensitivity + (sqrt survival_probability) ) )
     ; set immune_susceptibility abs ((immune_susceptibility - (sqrt chemo_sensitivity)))
      ]

      ;if pd1 = 1 [
      ;set immuno_sensitivity abs ( ( immuno_sensitivity - (sqrt chemo_sensitivity)) )
     ; set chemo_sensitivity abs ((chemo_sensitivity - (sqrt immuno_sensitivity)) )

;set doublebind_doublebind_immuno_sensitivity random-float 1
;set survival_probability rate_survival * 0.001
;set immune_susceptibility (1 - (doublebind_doublebind_immuno_sensitivity * (1 - survival_probability)))

;set CTLA4 random 2 ;commenting these out causes the population to crash
;set PD1 random 2


end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
647
448
-1
-1
13.0
1
10
1
1
1
0
1
1
1
-16
16
-16
16
1
1
1
ticks
30.0

BUTTON
42
81
167
114
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
46
40
112
73
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
654
317
826
350
cc
cc
0
10
2.0
1
1
NIL
HORIZONTAL

PLOT
858
10
1169
152
population
Generations
Population
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"Greens" 1.0 0 -14439633 true "" "if count turtles > 0 [plot (green-count / count turtles) * 100 ]"
"Blues" 1.0 0 -14070903 true "" "if count turtles > 0 [plot (blue-count / count turtles) * 100 ]"

BUTTON
83
193
175
226
single go
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
859
157
1092
394
PD1 CTLA4 +
NIL
NIL
0.0
10.0
0.0
100.0
true
true
"" ""
PENS
"pd1" 1.0 0 -16777216 true "" "if count turtles > 0 [ plot (pd1pop / count turtles ) * 100]"
"ctla4" 1.0 0 -1184463 true "" "if count turtles > 0 [plot  (CTLA4pop / count turtles ) * 100]"

PLOT
867
400
1242
574
PD1 CTLA4 -
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"ctla4_inactive" 1.0 0 -16777216 true "" "if count turtles > 0 [plot ctla4-pop / population]"
"pd1_inactive" 1.0 0 -7500403 true "" "if count turtles > 0 [plot pd1-pop / population]"

SLIDER
660
87
832
120
apoptosis_rate
apoptosis_rate
0
100
5.0
5
1
NIL
HORIZONTAL

MONITOR
210
549
290
594
population
population
17
1
11

SLIDER
658
201
833
234
doublebind_immune
doublebind_immune
0
1
1.0
.1
1
NIL
HORIZONTAL

SWITCH
662
240
806
273
chemo_switch
chemo_switch
1
1
-1000

SWITCH
661
278
806
311
immuno_switch
immuno_switch
1
1
-1000

MONITOR
1175
10
1232
55
Blue
blue-count
17
1
11

MONITOR
1177
60
1234
105
Green
green-count
17
1
11

MONITOR
35
360
177
405
Maximum Population
max_pop
17
1
11

PLOT
176
612
376
762
Both Inactive
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "if count turtles > 0 [plot  both_inactive / population]"

BUTTON
32
237
122
270
Prliferate
setupPop
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
205
509
378
542
NIL
ask turtles [proliferate]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
207
473
408
506
NIL
ask turtles [reach_capacity]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
660
124
832
157
doublebind_chemo
doublebind_chemo
0
1
1.0
.1
1
NIL
HORIZONTAL

SLIDER
661
163
833
196
chemo_dose
chemo_dose
0
1
1.0
.01
1
NIL
HORIZONTAL

SLIDER
659
48
831
81
immune_susceptibility
immune_susceptibility
0
1
0.26
.01
1
NIL
HORIZONTAL

SLIDER
658
10
830
43
max_pro
max_pro
0
1
0.65
.01
1
NIL
HORIZONTAL

TEXTBOX
9
505
176
631
Average proliferation/survival rate \nAs chemo dose increases, there is more selection for survival\nWhen immune susceptibility is high, the Pd1 and CTLA4 activated cells are selected for
11
0.0
1

TEXTBOX
57
133
207
161
Change plots to percentages
11
0.0
1

PLOT
428
640
628
790
plot 2
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"pd1" 1.0 0 -16777216 true "" "if population > 0 [plot  pd1pop]"
"ctla4" 1.0 0 -7500403 true "" "if population > 0 [plot  ctla4pop]"

PLOT
827
587
1244
840
plot 1
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"pd1 only" 1.0 0 -13791810 true "" "plot justpd1"
"ctla4 only" 1.0 0 -2139308 true "" "plot justctla4"
"Both active" 1.0 0 -987046 true "" "plot bothpop"

PLOT
427
476
627
626
Average Surv Rate 
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot (averagesurvivalrate * 10)"

MONITOR
35
304
181
349
Average Survival Rate
averagesurvivalrate * 100
5
1
11

SLIDER
657
357
829
390
chemo_start
chemo_start
0
200
50.0
1
1
NIL
HORIZONTAL

SLIDER
658
449
830
482
immuno_start
immuno_start
0
400
50.0
1
1
NIL
HORIZONTAL

SWITCH
655
497
801
530
timed_therapy
timed_therapy
1
1
-1000

SWITCH
656
544
834
577
population_therapy
population_therapy
0
1
-1000

SLIDER
657
404
829
437
population_start
population_start
0
1000
1000.0
50
1
NIL
HORIZONTAL

MONITOR
35
414
208
459
Average Proliferation Rate
averageproliferationrate
17
1
11

PLOT
627
591
827
741
Avg Proliferation Rate 
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot (averageproliferationrate * 10)"

@#$#@#$#@
## WHAT IS IT?

This model is a basic implementation of a double-bind using chemo and immuno-therapy. It assumes that cells forming resistance to one treatment directly affects their ability to form resistance against the other treatment. 

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

Generally we have had realistic results with the following parameters:
max_pro 0.65
immune_susceptibility 0.26
apoptosis_rate 5
cc (carrying capacity) : 2. If you're running the simulation on a computer more powerful than a MacBook Pro it may be appropriate to increase the carrying capacity. 
***
After that, you can change the values of the other parameters accordingly.
doublebind_chemo : the extent to which forming resistance to chemotherapy affects the cells' ability to form resistance to immunotherapy
doublebind_immune: inverse of doublebind_chemo
chemo_dose: the percentage strength of chemotherapy (100 = full strength, 50 = half strength, etc)
chemo/immuno_switch : administers individual therapies to the simulation. NOTE: timed therapy and population therapy will automatically switch the buttons on when appropriate, so have then set to "off" initially if using those therapy formats
chemo_start: the amount of ticks to begin chemotherapy
immuno_start: the amount of ticks to start immunotherapy
population_start: the minimum threshold for the population size to be in order to start treatment. Think of it as a way to allow the tumor to grow before treatment. 
timed_therapy/population_therapy: boolean values that will let the user decide if they want to base the therapy timing off of the population size or the general amount of ticks 
## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.3
@#$#@#$#@
@#$#@#$#@
1.0
    org.nlogo.sdm.gui.AggregateDrawing 1
        org.nlogo.sdm.gui.StockFigure "attributes" "attributes" 1 "FillColor" "Color" 225 225 182 137 175 60 40
            org.nlogo.sdm.gui.WrappedStock "popCTLA4" "100" 0
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2400"/>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="doublebind_chemo" first="0.1" step="0.1" last="0.9"/>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="20"/>
    </enumeratedValueSet>
    <steppedValueSet variable="doublebind_immune" first="0.1" step="0.1" last="0.9"/>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2400"/>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ChemoDose">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0.7"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immuno_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ChemoDose">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0.7"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ChemoDose">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0.7"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ChemoDose">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0.7"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment apop 10" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ChemoDose">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0.7"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="test" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2400"/>
    <metric>count turtles</metric>
    <steppedValueSet variable="doublebind_immune" first="0" step="0.2" last="1"/>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="chemo_dose" first="0" step="0.2" last="1"/>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="max_pro" first="0.3" step="0.1" last="1"/>
    <steppedValueSet variable="doublebind_chemo" first="0" step="0.2" last="1"/>
    <steppedValueSet variable="immune_susceptibility" first="0.1" step="0.1" last="1"/>
    <enumeratedValueSet variable="immuno_switch">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="1/26/17" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2400"/>
    <metric>count turtles</metric>
    <steppedValueSet variable="doublebind_immune" first="0" step="0.5" last="1"/>
    <steppedValueSet variable="immune_susceptibility" first="0.3" step="0.3" last="0.9"/>
    <enumeratedValueSet variable="max_pro">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="3"/>
    </enumeratedValueSet>
    <steppedValueSet variable="chemo_dose" first="0" step="0.3" last="0.9"/>
    <steppedValueSet variable="doublebind_chemo" first="0" step="0.5" last="1"/>
  </experiment>
  <experiment name="experiment1/27/17" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="0.35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="doublebind_chemo" first="0" step="0.2" last="1"/>
    <steppedValueSet variable="immune_susceptibility" first="0" step="0.2" last="0.1"/>
    <steppedValueSet variable="doublebind_immune" first="0" step="0.2" last="1"/>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="1/29/17" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2400"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <steppedValueSet variable="chemo_dose" first="0" step="0.2" last="1"/>
    <steppedValueSet variable="doublebind_immune" first="0" step="0.2" last="1"/>
    <steppedValueSet variable="doublebind_chemo" first="0" step="0.2" last="1"/>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="immune_susceptibility" first="0" step="0.5" last="1"/>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2400"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="30"/>
    </enumeratedValueSet>
    <steppedValueSet variable="chemo_dose" first="0" step="0.5" last="1"/>
    <steppedValueSet variable="doublebind_chemo" first="0" step="0.5" last="1"/>
    <steppedValueSet variable="doublebind_immune" first="0" step="0.5" last="1"/>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2/8/17" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2400"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0"/>
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0"/>
      <value value="0.3"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="0"/>
      <value value="0.3"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0"/>
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>count pd1</metric>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="3/16/17" repetitions="2" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1500"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <metric>count pd1pop</metric>
    <metric>count ctla4pop</metric>
    <metric>count averagesurvivalrate</metric>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.26"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ImmuneSusTest" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="immune_susceptibility" first="0.2" step="0.02" last="0.3"/>
  </experiment>
  <experiment name="ImmuneSusTestSecond" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.24"/>
      <value value="0.25"/>
      <value value="0.26"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ImmuneSusTestThird" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.26"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="NewData" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1500"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <metric>count turtles</metric>
    <metric>count pd1pop</metric>
    <metric>count ctla4pop</metric>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.26"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Final DATA" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1500"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.26"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="db 0" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1500"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="immuno_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.26"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="3/26 chemo test" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1500"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.26"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="immune - chemo no overlap" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles</metric>
    <metric>averagesurvivalrate</metric>
    <metric>averageproliferationrate</metric>
    <enumeratedValueSet variable="doublebind_chemo">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="apoptosis_rate">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="doublebind_immune">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_pro">
      <value value="0.65"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_dose">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cc">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immune_susceptibility">
      <value value="0.26"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_switch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuno_start">
      <value value="50"/>
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="chemo_start">
      <value value="25"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population_start">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timed_therapy">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population_therapy">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
