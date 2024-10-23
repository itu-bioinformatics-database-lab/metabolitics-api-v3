from marshmallow import Schema, fields
from flask_marshmallow import Marshmallow

from .app import app
from .models import User, Analysis

ma = Marshmallow(app)

#
class AnalysisInputSchema(Schema):
    study_name = fields.String(required=True)
    public = fields.Boolean(required=True)
    # concentration_changes = fields.Dict(required=True)
    analysis = fields.Dict(required=True)
    group = fields.String(required=True)
    isMapped = fields.Dict(required=False)
    disease = fields.Integer(required=True)
    metabolites = fields.List(required=False, cls_or_instance=fields.String())

class AnalysisInputSchema2(Schema):
    study_name = fields.String(required=True)
    public = fields.Boolean(required=True)
    # concentration_changes = fields.Dict(required=True)
    analysis = fields.Dict(required=True)
    group = fields.String(required=True)
    disease = fields.Integer(required=True)
    isMapped = fields.Dict(required=False)
    email = fields.String(required=True)


class PasswordChangeSchema(Schema):
    old_password = fields.String(required=True)
    new_password = fields.String(required=True)

#
class UserSchema(ma.ModelSchema):
    name = fields.String(required=True)
    surname = fields.String(required=True)
    email = fields.Email(required=True)
    affiliation = fields.String(required=True)
    password = fields.String(required=True, load_only=True)
    # confirmPassword = fields.String(required=True, load_only=True)

    class Meta:
        model = User
        exclude = ('analysis', )


class AnalysisSchema(ma.ModelSchema):
    results = fields.Dict()
    visualization = fields.Dict()

    # type = fields.String(required=False)

    class Meta:
        model = Analysis
        exclude = ('user', )


class PathwayChangesScheme(Schema):
    pathway = fields.String(required=True)
    change = fields.Integer(required=True)
    qualifier = fields.String(allow_none=True)
    amount = fields.Number(allow_none=True)
