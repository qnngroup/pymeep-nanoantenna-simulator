;;Calculates response of a triangular nanoparticle with a non-zero
;;radius of curvature at the tip apex.

;; Meep Unit Used: 1 um
;; Beam incident along z-direction.
;; Particle a 2-D triangle in xy with finite thickness.

;;Fixed Simulation Settings (Not able to be set by user)
(set! ensure-periodicity false) ;;Don't make shapes periodic (this has no effect on field periodicity)


;;EPS Averaging
(define-param avgeps? false)
(set! eps-averaging? avgeps?)

;;Structure Settings
(define-param structure? false) ; Sets the structure on or off for reflectivity calculation.
 
(define-param visualize? false) ; Visualization mode assists in viewing the drawing to make sure it is correct.  Shapes are filled with constant epsilon.
										;The epsilon values are arbitrary for color definition only.
(define-param sample 10) ; Sampling rate (Samples-per-period of highest frequency term in simulation)

;;Source Terms
(define-param fcen 3.0)
(define-param df 6.0)
(define-param amp 1)
(define-param nfreq 100) ; number of frequencies at which to compute flux


;;Simulation Space Definitions
(define-param dpml 0.1) ; PML thickness (applied below in Z direction only!)
						; Periodic in X and Y (see k-point setting in BC section below)
(define-param sx 0.52) ; size of x-dimension
(define-param sy 0.52) ; size of y-dimension
(define-param sz 1.0) ; size of z-dimension
(define-param res 300) ; Resolution: points per meep-unit
                       ; e.g. if using 1 meep unit = 1 um,
                       ; 1/res um spacing

;;Structure Settings
(define-param metal_thickness 0.02) ; Thickness of particle
(define-param contact_thickness 0.06) ; contact material thickness (e.g. ITO)
(define-param altitude .26) ; Altitude (height) of the triangle.
(define-param r_curvature 0.005);;Radius of curvature
(define-param contact_eps_r 3.06188) ;; real of epsilon for contact material
(define-param contact_sig_d 0.179556) ;; absorption coefficient of contact material
                                      ;;(meep-style)

;;Calculate actual size in Z of simulation space including the pmls!
(define sZ (+ (* 2.0 dpml) sz))

;;;;;;;;;;;;;;;;;;;;;;
;;Material Properties
;;;;;;;;;;;;;;;;;;;;;;

;;Gold - Adapted by Bala (see Bala's Blog)
(define myAu (make dielectric (epsilon 1)
(polarizations
 (make polarizability
(omega 1e-20) (gamma 0.042747) (sigma 4.0314e+41))
(make polarizability
(omega 0.33472) (gamma 0.19438) (sigma 11.363))
(make polarizability
(omega 0.66944) (gamma 0.27826) (sigma 1.1836))
(make polarizability
(omega 2.3947) (gamma 0.7017) (sigma 0.65677))
(make polarizability
(omega 3.4714) (gamma 2.0115) (sigma 2.6455))
(make polarizability
(omega 10.743) (gamma 1.7857) (sigma 2.0148))
)))
;Material used is Au from Rakic et al.,Applied Optics (1998)
;Plasma Angular Frequency (and plasma wave vector,kp) in normalized units=6.3493

(define mySapphire
      (make dielectric (epsilon 1.0)
            (E-susceptibilities 
             (make lorentzian-susceptibility
               (frequency (/ 1 (sqrt .0052799261))) (gamma 0.0) (sigma 1.43134930))
             (make lorentzian-susceptibility
               (frequency (/ 1 (sqrt 0.0142382647))) (gamma 0.0) (sigma 0.65054713))
             (make lorentzian-susceptibility
               (frequency (/ 1 (sqrt 325.017834))) (gamma 0.0) (sigma 5.3414021))
			 )))

;;ITO Defined as an absorbing thin-film.  
(define myITO (make dielectric
				(epsilon contact_eps_r)
				(D-conductivity contact_sig_d)
				)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Define object materials
;;;;;;;;;;;;;;;;;;;;;;;;;;
(define air)
(define substrate)
(define particle_material)
(define contact_material)

;; If visualization mode on, fill with placeholder epsilons.
;; If not, then fill with desired epsilons.
(if visualize?
	(begin ;; Visualization Mode
	 (set! air (make dielectric (epsilon 1)))
	 (set! substrate (make dielectric (epsilon 2)))
	 (set! particle_material (make dielectric (epsilon 3)))
	 (set! contact_material (make dielectric (epsilon 4)))
	 )

	(begin ;; Simulation Mode (Define actual Values)
	 
	 ;;Air, or vacuum:
	 (set! air (make dielectric (epsilon 1)))

	 ;;Substrate Material
	 (set! substrate (make dielectric (epsilon (expt 1.7522 2))))

	 ;;Particle Material
	 (set! particle_material myAu)

	 (set! contact_material myITO)
	 
	 )

	;;End If
	)


;;;;;;;;;;;;;;;;;;;;;;
;; Geometry
;;;;;;;;;;;;;;;;;;;;;;

;;Objective: create a single triangle with a rounded tip.  
(define base (* altitude 0.75))

(define theta (atan (/ 8 3))) ;;Angle of bottom right corner of triangle

(define r_top (* r_curvature (tan theta))) ;;length of side from apex
;;to corner where intersecting circle.

(define height_cutout (* r_top (cos (- (/ pi 2.0) theta))));;height of the cutout region of original
;;triangle

(define half_base_cutout (* r_top (sin (- (/ pi 2.0) theta))));;half of the base width
;;of cutout region on top of triangle

(define altitude_adjusted (- altitude height_cutout)) ;;New altitude of cut-off triangle

(define side (sqrt (+ (expt (- (* 0.5 base) half_base_cutout) 2.0) (expt altitude_adjusted 2.0))))
;;Side length of parallelograms to cutout each side of box to make triangles.

(define cutout_center_offset (+ half_base_cutout (/ (* (cos theta) side) 2.0))) ;;offset used to
;;determine exact position of center for parallelogram cutouts.

;;Set the lattice size of the entire simulation cell:
(set! geometry-lattice (make lattice (size sx sy sZ)))
;;If the structure is on, draw the substrate, contact material and nano-triangle.
;;If not, then just fill with air for field reference.
(set! geometry (if structure?
				  (list

				   ;;Bottom Triangle
				   (make block
				   	 (center 0
				   			 0
				   			 (- (/ metal_thickness -2.0) contact_thickness)
				   			 )
				   	 (size base altitude_adjusted metal_thickness)
				   	 (material particle_material)
				   	 )
				   (make block
				   	 (center (+ (/ base 4.0) cutout_center_offset)
				   			 0
				   			 (- (/ metal_thickness -2.0) contact_thickness)
				   			 )
				   	 (e2 (/ base -2.0) altitude 0)
				   	 (size (/ base 2.0) side metal_thickness)
				   	 (material air)
				   	 )
				   (make block
				   	 (center (- (/ base -4.0) cutout_center_offset)
				   			 0
				   			 (- (/ metal_thickness -2.0) contact_thickness)
				   			 )
				   	 (e2 (/ base 2.0) altitude 0)
				   	 (size (/ base 2.0) side metal_thickness)
				   	 (material air)
				   	 )

				   ;;Round Bottom Triangle
	               ;; -- Note one has to be careful of y-positioning to ensure it
                   ;;    is properly aligned with blunted triangle we have formed.
				   (make cylinder
				   	 (center 0
				   			 (- (* altitude_adjusted 0.5) (* (cos theta) r_curvature))
				   			 (- (/ metal_thickness -2.0) contact_thickness))
				   	 (radius r_curvature)
				   	 (height metal_thickness)
				   	 (material particle_material)					
				   	 )

				   ;;Contact Layer
				   (make block (center 0 0 (/ contact_thickness -2.0))
						 (size infinity infinity contact_thickness)
						 (material contact_material)
						 )

				  ;;Substrate				   
				  (make block (center 0 0 (* sZ 0.25)) 
                    (size infinity infinity (* sZ 0.5))
                    (material substrate))
				   
              )
              (list (make block (center 0 0 0) (size infinity infinity infinity)
						  (material air))
				   ;;Include this if you want the substrate included in the extinction calculation.
			       ;; (make block (center 0 0 (* sZ 0.25)) 
                   ;;  (size infinity infinity (* sZ 0.5))
					;;  (material substrate))
			  )
			)
)

;;;;;;;;;;;;;;;;;;;;;;
;; PML Layers
;;;;;;;;;;;;;;;;;;;;;;

;;Define the PML layers now that the epsilon values are all set
(set! pml-layers (list (make pml (direction Z) (thickness dpml))))

;;;;;;;;;;;;;;;;;;;;;;
;; Symmetry
;;;;;;;;;;;;;;;;;;;;;;

;;Exploit symmetry (if folded along the y-axis, so direction is X)
(set! symmetries (list (make mirror-sym (direction X))
					   )
	  )

;;;;;;;;;;;;;;;;;;;;;;;;
;; Boundary Conditions
;;;;;;;;;;;;;;;;;;;;;;;;
(set! k-point (vector3 0 0 0)) ;;kvector of zero is used to
                               ;;keep same phase in all the units cells.


;;;;;;;;;;;;;;;;
;; Sources
;;;;;;;;;;;;;;;;
(set! resolution res)

;Source definition
(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component Ey)
				 (center 0 0 (* sz -0.125))
				 (size sx sy 0))))

;;;;;;;;;;;;;;;;
;; Monitors
;;;;;;;;;;;;;;;;

;;Transmission Flux Monitor
(define trans                                                
      (add-flux fcen df nfreq
                    (make flux-region
					  (center 0 0 (* sz 0.25))
					  (size sx sy 0)
					  (direction Z))))

;Reflected Flux Monitor
(define refl ; reflected flux
      (add-flux fcen df nfreq
                 (make flux-region
                   (center 0 0 (* sz -0.25))
				   (size sx sy 0))))

;;Field Monitor
(define monitor-xy
  (volume 
   (center 0 0 (- (* metal_thickness -0.5) contact_thickness))
   (size sx sy 0)
   )
  )

;;Diffraction plane monitor
(define diff-xy
  (volume 
   (center 0 0 (* sz 0.5))
   (size sx sy 0)
   )
  )
 
;; If the structure is there, then load the reference
;; flux from the no structure case.  This ensures
;; a correct calculation of reflected flux.
(if (not visualize?)
	(if structure? (load-minus-flux "refl-flux" refl))		   
	)

;;;;;;;;;;;;;;;;;;;;
;; Run the sources
;;;;;;;;;;;;;;;;;;;;

;;If visualize is on, just run to see the epsilon configuration...
;;else, then run the full simulation!x
(if visualize? (run-until 0 (to-appended "t_vis" (at-beginning output-epsilon)))

	;;If structure, then record Ex, Ey and Ez
	(if structure? 
		(run-sources+  (stop-when-fields-decayed 10 Ey
												 (vector3 0 0 (* sz 0.25))
												 1e-3)

					   ;; For reference, take snapshot of eps.
					   (at-beginning output-epsilon)

					   ;; Sample Ex
					   (in-volume monitor-xy
								  (to-appended "tEx_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-x)))
					   ;;Sample Ey 
					   (in-volume monitor-xy
								  (to-appended "tEy_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-y)))
					   ;;Sample Ez
					   (in-volume monitor-xy
								  (to-appended "tEz_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-z)))

					   ;;Sample Diff Plane Ex
					   (in-volume diff-xy
								  (to-appended "tEx_diff_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-x)))

					   ;;Sample Diff Plane Ey
					   (in-volume diff-xy
								  (to-appended "tEy_diff_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-y)))

					   ;;Sample Diff Plane Ez
					   (in-volume diff-xy
								  (to-appended "tEz_diff_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-z)))

					   )

		;;If no structure, just record Ey.
		(run-sources+  (stop-when-fields-decayed 5 Ey
												 (vector3 0 0 (* sz 0.25))
												 1e-3)

					   ;;Just Ey when there is no nano-particle.  
					   (in-volume monitor-xy
								  (to-appended "bEy_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-y)))

					   ;;Background Diff Plane Ey
					   (in-volume diff-xy
								  (to-appended "bEy_diff_xy"
											   (at-every
												(/ 1 (+ fcen (* 0.5 df)) sample)
												output-efield-y)))
			)		   
					
		)

	)
                            

;;If not just visualizing, then save the reflected and transmitted data.
(if (not visualize?)
	(begin
	 ;; Save reflected flux if there is no nano
	 ;; particle present.  This is for reference when
	 ;; calculating correct reflected flux with the particle.
	 (if (not structure?)
		 (save-flux "refl-flux" refl))
	 
	 ;;Now output the fluxes to grep
	 (display-fluxes refl trans)

	 )
	)





