{-# LANGUAGE PackageImports #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DataKinds #-}

import Data.Bits
import Data.Sized.Matrix as M
import Data.Sized.Sampled as SMP
import Data.Sized.Signed as S
import Data.Ix

type instance Index (Signed ix)  = Int

instance forall ix . (Size ix) => Ix (Signed ix) where
    range     (l, u)    = [l .. u]
    inRange   (l, u) v  =  (l <= v) && (v <= u)
    index     (l, u) v | inRange (l,u) v = fromIntegral (v - l)
                       | otherwise       = error "Error in Ix array index"
    rangeSize (l, u)   | l <= u           = fromIntegral $ (toInteger u) - (toInteger l) + 1
                       | otherwise       = 0

instance forall ix . (Size ix) => Size (Signed ix) where
    size         = const s
       where s  = fromIntegral $ toInteger (maxBound :: Signed ix)
    addIndex v n =  v + (fromIntegral n)  -- fix bounds issues
    toIndex v    = fromIntegral v


type F16 = Sampled S16 X32 -- Fixed point representation of 32 bit signed number; 15 bits decimal + 16 bit fractional (15.16)

rad :: [Double] -- atan(2 ^ (-i)) angles in radians
rad = [ atan ( 2 ^^ (-i) ) | i <- [ 0 .. ] ]

radF :: [F16] -- convert to fixed point representation F16
radF = map (mkSampled . toRational) $ rad

shiftR' :: F16 -> Int -> F16 -- convert to signed number, shift, then convert back to sampled
shiftR' n i = SMP.fromMatrix (S.toMatrix (shiftR signed_n i))
 where signed_n = S.fromMatrix (SMP.toMatrix n)

scale :: Int -> Double
scale i = 1 / sqrt ( 1 +  2 ^^ (( -2 ) * i ))

sf :: [Double] -- scaling factors
sf = list_sf 1 ( scale 0 )
  where list_sf :: Int -> Double -> [Double]
        list_sf i n = n : list_sf ( i + 1 ) ( scale i * n )

sfF :: [F16]  -- convert to fixed point representation F16
sfF = map (mkSampled . toRational) $ sf

calc_xy :: Int -> F16 -> F16 -> F16 -> Matrix X2 F16 -- next x and y
calc_xy i sign x y = matrix (x' : y' : [])
 where x' = x - sign * (shiftR' y i )
       y' = y + sign * (shiftR' x i )

iter_vec :: (Int, F16, Matrix X2 F16) -> F16 -> (Int, F16, Matrix X2 F16) -- vectoring mode
iter_vec (i, z, xy) alpha
 | y > 0 = ( i+1 , z + alpha, calc_xy i (-1) x y)
 | otherwise = ( i+1 , z - alpha, calc_xy i 1 x y)
 where x = xy ! 0
       y = xy ! 1

iter :: (Int, F16, Matrix X2 F16) -> F16 -> (Int, F16, Matrix X2 F16) -- rotation mode
iter (i, z, xy) alpha
  | z < 0 = ( i+1 , z + alpha, calc_xy i (-1) x y)
  | otherwise = ( i+1 , z - alpha, calc_xy i 1 x y)
  where x = xy ! 0
        y = xy ! 1

cordic :: F16 -> Int -> Matrix X2 F16 -- input is in radian, return matrix of cos(input) and sin(input) for n_iter of iterations
cordic input n_iter = matrix ( c_s ! 0 * sf : c_s ! 1 * sf : [])
 where init = ( 0, input, matrix (1 : 0 : []) ) -- initialize to x0 = 1, y0 = 0; to avoid the last multiplication by scaling factor initialize x to K = 0.607252935
       (i, z, c_s) = foldl iter init $ take n_iter radF
       sf = sfF !! i

cordic_vec :: F16 -> F16 -> Int -> F16 -- input is in radian, return matrix of cos(input) and sin(input) for n_iter of iterations
cordic_vec x y n_iter = z
 where init = ( 0, 0, matrix (x : y : []) ) -- initialize to x0 = 1, y0 = 0; to avoid the last multiplication by scaling factor initialize x to K = 0.607252935
       (i, z, c_s) = foldl iter_vec init $ take n_iter radF
       --sf = sfF !! i
