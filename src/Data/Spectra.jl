import CSV, DataFrames

kopinoledpanel() = MeasuredSpectrum(CSV.read(joinpath(@__DIR__, "OLED Spectrum Kopin panel.csv")))
export kopinoledpanel
