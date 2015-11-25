/*
 *  Models & Collections
 */
Relation = Backbone.RelationalModel.extend({
  defaults: {
    id: '',
    text: '',
    selected: false
  },

  relations: [{
    type: 'HasMany',
    key: 'children',

    relatedModel: 'Relation',
    collectionType: 'RelationList',

    reverseRelation : {
      key : 'parentRelation',
      includeInJSON: false,
    }
  }],

  get_selected: function() {
    /* Try to find an instance of 'this' model type in the store */
    var model = Backbone.Relational.store.find( this, {"selected": true});
    console.log('model:', model);

    if ( !model && _.isObject( attributes ) ) {
      var coll = Backbone.Relational.store.getCollection( this );

      model = coll.find( function( m ) {
        return m.selected === true;
      });
    }
    return model;
  }

});

RelationList = Backbone.Collection.extend({
  model: Relation,
  url: '/api/v1/words',
});


/*
 * Views
 */
RelationView = Backbone.Marionette.ItemView.extend({
  template: _.template('<%= text %>'),
  tagName: 'a',
  className: 'list-group-item',

  events : {
    'mousedown' : 'mousedown',
  },

  /*
  findPersonInPeopleCollection: function(nameWeAreLookingFor) {
    var model = this.model;

    function findChildren(obj) {
      if (!obj.children) return [obj];
      var children = _.map(obj.children, function(child) {
        foundPeople.push( findChildren(child) );
      });
      return _.flatten(children);
    }

    var allPeople = _.flatten( model.collection.map(findChildren) );

    return _(allPeople).find(function(obj) {
      return obj.get('selected') == true;
    }
  },
  */

  mousedown : function(evt) {
    /*
     * 1. Set the current choice to the ID
     * 2. Set the backbutton reference if available
     * 3. Set the list updated to any children
     *
     */
    var children = this.model.get('children');
    Tree['convoChannel'].trigger('click', {'collection': children, 'choice': this.model});
  },
});


RelationCompositeView = Backbone.Marionette.CompositeView.extend({
  template: '#tree-template',
  templateHelpers: function() {
    return this.options.concepts
  },

  childView  : RelationView,
  childViewContainer: "ul",
  /*
  tagName   : 'div',
  className : 'paragraph',
  */

  ui: {
    'c1': '#c1',
    'relation': '#relation',
    'c2': '#c2',
    'c1_not_correct': '#c1_not_correct',
    'c1_remove_icon': '#c1_remove_icon',
    'c2_not_correct': '#c2_not_correct',
    'c2_remove_icon': '#c2_remove_icon',
  },

  events : {
    'mousedown @ui.relation': 'resetRelationship'
  },

  resetRelationship: function(evt) {
    Tree['convoChannel'].trigger('back', this.options);
  },

  onRender : function() {
    var self = this;
    var choice = this.options['choice'];
    var concepts = this.options['concepts'];
    /*
     * console.log('[RelationCompositeView onRender] Choice:', choice);
     */

    if(choice) {
      this.ui.relation.removeClass('disabled').text( choice.get('text') );
    }

    if(choice || this.collection.parentRelation) {
      this.ui.relation.removeClass('disabled');
      self.ui.relation.addClass('relation-go-back');
    }

    function get_stype_word(stype) {
      if (stype === 'g')
        stype_word = "gene";
      else if (stype === "d")
        stype_word = "disease";
      else if (stype === "c")
        stype_word = "drug";
    return stype_word;
    };

    var c1TimeoutId;
    c1_stype_word = get_stype_word(concepts['c1'].type);
    c2_stype_word = get_stype_word(concepts['c2'].type);

    // this.ui.c1.fadeIn(700); //TODO I want this but not every time the menu is clicked.
    this.ui.c1_remove_icon.addClass("fa fa-times-circle").css({'font-size': '25px', 'color': 'grey'});
    this.ui.c1_not_correct.css({"color":color_find(concepts['c1'].type)});
    this.ui.c1_not_correct.html('<h3>is not a '+ c1_stype_word +'.</h3>')
    this.ui.c1.css({"background-color":color_find(concepts['c1'].type),"min-height":"220px"});
    this.ui.c1.hover(function() {
      if (!c1TimeoutId) {
        c1TimeoutId = window.setTimeout(function() {
          c1TimeoutId = null;
          // self.ui['c1'].addClass('not-correct-concept');
          self.ui['c1_remove_icon'].css("color", "red")
          self.ui['c1'].css("background-color", "red");
          self.ui['c1_not_correct'].html('<h3>is not a '+ c1_stype_word +'.</h3>').css("color", "black");

        }, 600);
      }
    }, function() {
      if (c1TimeoutId) {
        window.clearTimeout(c1TimeoutId);
        c1TimeoutId = null;
        self.render();
      } else {
        console.log('elsed');
        self.render();
      }
    });

    function color_find(relationship_type) {
      var color;
      if (relationship_type === "g") { color = "#B1FFA8"; };
      if (relationship_type === "d") { color = "#d1f3ff"; };
      if (relationship_type === "c") { color = "#ffd1dc"; };
      return color;
    };

    var c2TimeoutId;
    // this.ui.c2.fadeIn(700);
    this.ui.c2_remove_icon.addClass("fa fa-times-circle").css({'font-size': '25px', 'color': 'grey'});
    this.ui.c2_not_correct.css({"color":color_find(concepts['c2'].type)});
    this.ui.c2_not_correct.html('<h3>is not a '+ c2_stype_word +'.</h3>')
    this.ui.c2.css({"background-color":color_find(concepts['c2'].type),"min-height":"220px"});
    this.ui.c2.hover(function() {
      stype_word = get_stype_word(concepts['c2'].type);
      if (!c2TimeoutId) {
        c2TimeoutId = window.setTimeout(function() {
          c2TimeoutId = null;
          // self.ui['c1'].addClass('not-correct-concept');
          self.ui['c2_remove_icon'].css("color", "red")
          self.ui['c2'].css("background-color", "red");
          self.ui['c2_not_correct'].html('<h3>is not a '+ c2_stype_word +'.</h3>').css("color", "black");
        }, 600);
      }
    }, function() {
      if (c2TimeoutId) {
        window.clearTimeout(c2TimeoutId);
        c2TimeoutId = null;
        self.render();
      } else {
        console.log('elsed');
        self.render();
      }
    });

  }
});

Tree = new Backbone.Marionette.Application();
