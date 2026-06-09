"""Schemas for tabular data used in the workflow."""

from pandera import Field, check
from pandera.pandas import DataFrameModel
from pandera.typing.geopandas import GeoSeries
from pandera.typing.pandas import Index, Series
from shapely.geometry import Point


class PointLayout(DataFrameModel):
    class Config:
        coerce = True
        strict = True

    id: Index[int] = Field(unique=True)
    "Unique ID for this layout point."
    techs: Series[str] = Field()
    "Technology type"
    lat: Series[float] = Field()
    "Latitude"
    lon: Series[float] = Field()
    "Longitude"
    capacity: Series[float] = Field()
    "Installed capacity"


class Shapes(DataFrameModel):
    class Config:
        coerce = True
        strict = False

    shape_id: Series[str] = Field(unique=True)
    "Unique ID for this shape."
    country_id: Series[str] = Field()
    "ISO alpha-3 code."
    shape_class: Series[str] = Field(isin=["land", "maritime"])
    "Shape classifier"
    geometry: GeoSeries[Point] = Field()
    "Shape polygon."

    @check("geometry")
    def geom_not_empty(cls, geom):
        return geom.notna().all() & (~geom.is_empty).all() & geom.is_valid.all()
