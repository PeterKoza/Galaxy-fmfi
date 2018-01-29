define("mvc/collection/collection-model",["exports","mvc/dataset/dataset-model","mvc/base-mvc","utils/localization"],function(t,e,n,i){"use strict";function o(t){return t&&t.__esModule?t:{default:t}}Object.defineProperty(t,"__esModule",{value:!0});var s=o(e),l=o(n),r=(o(i),{defaults:{model_class:"DatasetCollectionElement",element_identifier:null,element_index:null,element_type:null},_mergeObject:function(t){return _.extend(t,t.object,{element_id:t.id}),delete t.object,t},constructor:function(t,e){t=this._mergeObject(t),this.idAttribute="element_id",Backbone.Model.apply(this,arguments)},parse:function(t,e){var n=t;return n=this._mergeObject(n)}}),c=Backbone.Model.extend(l.default.LoggableMixin).extend(r).extend({_logNamespace:"collections"}),a=Backbone.Collection.extend(l.default.LoggableMixin).extend({_logNamespace:"collections",model:c,toString:function(){return["DatasetCollectionElementCollection(",this.length,")"].join("")}}),d=s.default.DatasetAssociation.extend(l.default.mixin(r,{url:function(){return this.has("history_id")?Galaxy.root+"api/histories/"+this.get("history_id")+"/contents/"+this.get("id"):(console.warn("no endpoint for non-hdas within a collection yet"),Galaxy.root+"api/datasets")},defaults:_.extend({},s.default.DatasetAssociation.prototype.defaults,r.defaults),_downloadQueryParameters:function(){return"?to_ext="+this.get("file_ext")+"&hdca_id="+this.get("parent_hdca_id")+"&element_identifier="+this.get("element_identifier")},constructor:function(t,e){this.debug("\t DatasetDCE.constructor:",t,e),r.constructor.call(this,t,e)},hasDetails:function(){return this.elements&&this.elements.length},toString:function(){return["DatasetDCE(",this.get("element_identifier"),")"].join("")}})),u=a.extend({model:d,toString:function(){return["DatasetDCECollection(",this.length,")"].join("")}}),h=Backbone.Model.extend(l.default.LoggableMixin).extend(l.default.SearchableModelMixin).extend({_logNamespace:"collections",defaults:{collection_type:null,deleted:!1},collectionClass:a,initialize:function(t,e){this.debug(this+"(DatasetCollection).initialize:",t,e,this),this.elements=this._createElementsModel(),this.on("change:elements",function(){this.log("change:elements"),this.elements=this._createElementsModel()})},_createElementsModel:function(){this.debug(this+"._createElementsModel",this.collectionClass,this.get("elements"),this.elements);var t=this.get("elements")||[];this.unset("elements",{silent:!0});var e=this;return _.each(t,function(t,n){_.extend(t,{parent_hdca_id:e.get("id")})}),this.elements=new this.collectionClass(t),this.elements},toJSON:function(){var t=Backbone.Model.prototype.toJSON.call(this);return this.elements&&(t.elements=this.elements.toJSON()),t},inReadyState:function(){var t=this.get("populated");return this.isDeletedOrPurged()||t},hasDetails:function(){return 0!==this.elements.length},getVisibleContents:function(t){return this.elements},parse:function(t,e){var n=Backbone.Model.prototype.parse.call(this,t,e);return n.create_time&&(n.create_time=new Date(n.create_time)),n.update_time&&(n.update_time=new Date(n.update_time)),n},delete:function(t){return this.get("deleted")?jQuery.when():this.save({deleted:!0},t)},undelete:function(t){return!this.get("deleted")||this.get("purged")?jQuery.when():this.save({deleted:!1},t)},isDeletedOrPurged:function(){return this.get("deleted")||this.get("purged")},searchAttributes:["name","tags"],toString:function(){return"DatasetCollection("+[this.get("id"),this.get("name")||this.get("element_identifier")].join(",")+")"}}),g=h.extend({collectionClass:u,toString:function(){return"List"+h.prototype.toString.call(this)}}),f=g.extend({toString:function(){return"Pair"+h.prototype.toString.call(this)}}),m=h.extend(l.default.mixin(r,{constructor:function(t,e){this.debug("\t NestedDCDCE.constructor:",t,e),r.constructor.call(this,t,e)},toString:function(){return["NestedDCDCE(",this.object?""+this.object:this.get("element_identifier"),")"].join("")}})),C=a.extend({model:m,toString:function(){return["NestedDCDCECollection(",this.length,")"].join("")}}),D=f.extend(l.default.mixin(r,{constructor:function(t,e){this.debug("\t NestedPairDCDCE.constructor:",t,e),r.constructor.call(this,t,e)},toString:function(){return["NestedPairDCDCE(",this.object?""+this.object:this.get("element_identifier"),")"].join("")}})),x=C.extend({model:D,toString:function(){return["NestedPairDCDCECollection(",this.length,")"].join("")}}),p=h.extend({collectionClass:x,toString:function(){return["ListPairedDatasetCollection(",this.get("name"),")"].join("")}}),b=g.extend(l.default.mixin(r,{constructor:function(t,e){this.debug("\t NestedListDCDCE.constructor:",t,e),r.constructor.call(this,t,e)},toString:function(){return["NestedListDCDCE(",this.object?""+this.object:this.get("element_identifier"),")"].join("")}})),j=C.extend({model:b,toString:function(){return["NestedListDCDCECollection(",this.length,")"].join("")}}),S=h.extend({collectionClass:j,toString:function(){return["ListOfListsDatasetCollection(",this.get("name"),")"].join("")}});t.default={ListDatasetCollection:g,PairDatasetCollection:f,ListPairedDatasetCollection:p,ListOfListsDatasetCollection:S}});
//# sourceMappingURL=../../../maps/mvc/collection/collection-model.js.map