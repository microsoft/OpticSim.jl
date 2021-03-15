const TEMP_REF = 20.0
const PRESSURE_REF = 1.0
const TEMP_REF_UNITFUL = TEMP_REF * u"°C"
export TEMP_REF, PRESSURE_REF

const DISPFORM_NAMES = ["Schott", "Sellmeier1", "Herzberger", "Sellmeier2", "Conrady", "Sellmeier3", "HandbookOfOptics1", "HandbookOfOptics2", "Sellmeier4", "Extended1", "Sellmeier5", "Extended2", "Extended3"]
const STATUS = ["Standard", "Preferred", "Obsolete", "Special", "Melt"]

const MIL_GLASSES = Dict{Int,Glass}()
const MODEL_GLASSES = Vector{Glass}(undef, 0)

"""
Special glass to represent air. Refractive index is defined to always be 1.0 for any temperature and pressure (other indices are relative to this).
"""
const Air = AirType()
