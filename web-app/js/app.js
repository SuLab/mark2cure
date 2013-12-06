var firefox = navigator.userAgent.toLowerCase().indexOf("firefox") > -1,
    explorer = navigator.userAgent.toLowerCase().indexOf("explorer") > -1;

if(firefox || explorer) {
  alert("Google Chrome or Safari are recommended.");
}

_.intersectionObjects = _.intersect = function(array) {
  var slice = Array.prototype.slice; // added this line as a utility
  var rest = slice.call(arguments, 1);
    return _.filter(_.uniq(array), function(item) {
      return _.every(rest, function(other) {
        //return _.indexOf(other, item) >= 0;
        return _.any(other, function(element) { return _.isEqual(element, item); });
      });
    });
};

vent.on('navigate:document:relationship', function(param) {
  var id = Number(param.doc_id),
      found_model = opts.collection.findWhere({id: id});

  if(found_model) {
    opts['model'] = found_model;
    app.main.show( new RelationGame(opts) );

  } else {
    //-- If the Doc wasn't already available
    var single_doc = new Document({id: id});

    single_doc.fetch({success : function() {
      opts['model'] = single_doc;

      app.main.show( new RelationGame(opts) );
    }, error : function() {
      vent.trigger('library', {});
    }});

  }
});

vent.on('navigate:document', function(param) {
  var id = Number(param.doc_id),
      found_model = opts.collection.findWhere({id: id});

  if(found_model) {
    opts['model'] = found_model;
    app.main.show( new Game(opts) );

  } else {
    //-- If the Doc wasn't already available
    var single_doc = new Document({id: id});

    single_doc.fetch({success : function() {
      opts['model'] = single_doc;

      if(opts.user.get('mturk') || window.aws && window.aws.assignment_id) {
        app.header.close();
        app.header.$el.hide();
        $('body').css({'padding-top':'20px'});
        app.footer.close();
      }

      app.main.show( new Game(opts) );
    }, error : function() {
      vent.trigger('library', {});
    }});

  }
});

vent.on('navigate:library', function() {
  //-- Do user checking here (prevents AMT from getting popup)
  if( !opts.user.authenticated() ) {
    app.header.close();
    app.main.close();
    app.footer.close();

    //-- This set the modal so it can't be closed
    opts.static = true;
    app.modal.show(   new FirstRun(opts) )
  } else {
    opts.collection.fetch();
  }

  //-- If first time / do training
  // var doc = this.collection.at(0)
  // Backbone.history.navigate( '#/document/'+ doc.id );
  // vent.trigger('navigate:document', doc.id );

  app.main.show( new Library(opts) );
});

vent.on('navigate:analytics', function(obj) {
  if(!obj['toggle'] || $('#network').length ) {
    app.analytics.close();
  } else {
    app.analytics.show( new Network(obj) );
  }
});
