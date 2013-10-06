define(['marionette', 'vent',
        //-- Models
        'models/Document', 'models/User',

        //-- Collections
        'collections/DocumentList',

        //-- Views
        'views/Header/Index', 'views/Footer/Index',
        'views/Main/Game/Index', 'views/Main/Library/Index',
        'views/Analytics/Network',
        'views/modals/FirstRun', 'views/modals/Complete', 'views/modals/Message', 'views/modals/Survey',
        'views/modals/Settings', 'views/modals/Instructions',

        //-- Utils
        'extensions/ModalRegion', 'extensions/FadeRegion', 'extensions/SlideRegion'],

   function(marionette, vent,
            //-- Models
            Document, User,

            //-- Collections
            DocumentList,

            //-- Views
            Header, Footer,
            Game, Library,
            Network,
            FirstRun, Complete, Message, Survey,
            Settings, Instructions,

            //--Utils
            ModalRegion, FadeRegion, SlideRegion) {
    'use strict';

    var app = new marionette.Application(),
        opts = {
          collection : new DocumentList(),
          user : User
        };

    app.addRegions({
      header    : 'header',
      analytics : SlideRegion.extend({el: "#analytics"}),
      main      : FadeRegion.extend({el: "#content"}),
      footer    : '#footer',
      modal     : ModalRegion,
    });

    app.addInitializer(function() {
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

      //-- If the user changes, refetch their docs
      this.listenTo(opts.user, 'change:api_key', opts.collection.fetch());

      app.header.show(  new Header(opts)  );
      app.main.show(    new Library(opts) );
      app.footer.show(  new Footer(opts)  );
    });

    vent.on('navigate:document', function(param) {
      var id = Number(param.doc_id),
          found_model = opts.collection.findWhere({id: id});

      if(found_model) {
        opts['model'] = found_model;
        app.main.show( new Game(opts) );

      } else {
        //-- If the Doc wasn't already available
        // var single_doc = new Document({id: id});

        // single_doc.fetch({success : function() {
        //   opts['model'] = single_doc;
        //   app.main.show( new Game(opts) );
        // }, error : function() {
        //   vent.trigger('library', {});
        // }});

      }
    });

    vent.on('navigate:library', function() {
      //-- Update URL
      Backbone.history.navigate( '#/library' );

      //-- Do user checking here (prevents AMT from getting popup)
      if( !opts.user.authenticated() ) {
        opts.user.save();

        app.header.close();
        app.main.close();
        app.footer.close();

        //-- This set the modal so it can't be closed
        opts.static = true;
        app.modal.show(   new FirstRun(opts) )
      }

      app.main.show( new Library(opts) );
    });

    vent.on('navigate:analytics', function(obj) {
      if(!obj['toggle'] || $('#network').length ) {
        app.analytics.close();
      } else {
        app.analytics.show( new Network(obj) );
      }
    });

    // vent.on('navigate:quest:search', function() {
    //   var quest = obj['quest']
    //   if(quest && quest.trim().length) {
    //     $.ajax({
    //       async: false,
    //       url: '/api/v1/quest/'+quest.trim(),
    //       success: function(jsonData) {
    //         var docs = []
    //         _.each(jsonData.objects, function(v) {
    //           var doc = new Document({id: v.document_id});
    //           doc.fetch();
    //           docs.push(doc);
    //         });
    //         viewOptions.collection = new DocumentList(docs);
    //         app.main.show( new Library(viewOptions) );
    //         if(docs.length==0) {
    //           vent.trigger('library', {});
    //         }
    //       },
    //       error: function() {
    //         Backbone.history.navigate( '#/library' );
    //         vent.trigger('library', {});
    //       }
    //     });
    // });

    //
    //-- Display Modals
    //
    vent.on('navigate:message', function() {
      app.modal.show( new Message(opts) );
    });

    vent.on('navigate:survey', function() {
      app.modal.show( new Survey(opts) );
    });

    vent.on('navigate:settings', function() {
      app.modal.show(   new Settings(opts) )
    });

    vent.on('navigate:instructions', function() {
      app.modal.show(   new Instructions(opts) )
    });

    //-- Oneoff modals really...
    vent.on('modal:show_complete', function() {
      viewOptions['swap'] = true;
      app.modal.show( new Complete(opts) );
    });

    vent.on('modal:welcome', function() {
      app.modal.show( new FirstRun(opts) );
    });

    vent.on('modal:close', function() {
      Backbone.history.navigate('#');
      $('#modal').modal('hide');
    });

    return app;

  }
);
