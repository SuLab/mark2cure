/**
 * Backbone-tastypie.js 0.2
 * (c) 2011 Paul Uithol
 *
 * Backbone-tastypie may be freely distributed under the MIT license.
 * Add or override Backbone.js functionality, for compatibility with django-tastypie.
 * Depends on Backbone (and thus on Underscore as well): https://github.com/documentcloud/backbone.
 */
(function( undefined ) {
	"use strict";

	// Backbone and underscore noConflict support. Save local reference to _ and Backbone objects.
	var _, Backbone;
	// CommonJS shim
	if ( typeof window === 'undefined' ) {
		_ = require( 'underscore' );
		Backbone = require( 'backbone' );
	}
	else {
		_ = window._;
		Backbone = window.Backbone;
	}

	Backbone.Model.prototype.url = function() {
		// Use the id if possible
		var url = this.id;

		// If there's no idAttribute, use the 'urlRoot'. Fallback to try to have the collection construct a url.
		// Explicitly add the 'id' attribute if the model has one.
		if ( !url ) {
			url = _.isFunction( this.urlRoot ) ? this.urlRoot() : this.urlRoot;
			url = url || this.collection && ( _.isFunction( this.collection.url ) ? this.collection.url() : this.collection.url );

			if ( url && this.has( 'id' ) ) {
				url = addSlash( url ) + this.get( 'id' );
			}
		}

		url = url && addSlash( url );

    console.log('protype url :: ', this, url);
		return url || null;
	};

	/**
	 * Return the first entry in 'data.objects' if it exists and is an array, or else just plain 'data'.
	 */
	Backbone.Model.prototype.parse = function( data ) {
		return data && data.objects && ( _.isArray( data.objects ) ? data.objects[ 0 ] : data.objects ) || data;
	};

	/**
	 * Return 'data.objects' if it exists.
	 * If present, the 'data.meta' object is assigned to the 'collection.meta' var.
	 */
	Backbone.Collection.prototype.parse = function( data ) {
		if ( data && data.meta ) {
			this.meta = data.meta;
		}

		return data && data.objects || data;
	};

	Backbone.Collection.prototype.url = function( models ) {

    var url = _.isFunction( this.urlRoot ) ? this.urlRoot() : this.urlRoot;
		// If the collection doesn't specify an url, try to obtain one from a model in the collection
		if ( !url ) {
			var model = models && models.length && models[ 0 ];
			url = model && ( _.isFunction( model.urlRoot ) ? model.urlRoot() : model.urlRoot );
		}
		url = url && addSlash( url );

		// Build a url to retrieve a set of models. This assume the last part of each model's idAttribute
		// (set to 'resource_uri') contains the model's id.
		if ( models && models.length ) {
			var ids = _.map( models, function( model ) {
				var parts = _.compact( model.id.split( '/' ) );
				return parts[ parts.length - 1 ];
			});
			url += 'set/' + ids.join( ';' ) + '/';
		}

		return url || null;
	};

	var addSlash = function( str ) {
		return str + ( ( str.length > 0 && str.charAt( str.length - 1 ) === '/' ) ? '' : '/' );
	};
})();
